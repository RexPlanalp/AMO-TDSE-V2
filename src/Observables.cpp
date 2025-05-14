#include "Observables.h"



void Observables::computeDistribution(int rank,const Basis& Basis, const TDSE& tdse,const TISE& tise, const Angular& angular)
{
    if (rank != 0)
    {
        return;
    }

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,Basis.getNbasis(),Basis.getNbasis(),2*Basis.getDegree() + 1};
    RadialMatrix::populateRadialMatrix(RadialMatrixType::S,S,Basis,false);

    PetscHDF5 viewer{PETSC_COMM_SELF,tdse.getOutputPath(),FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,Basis.getNbasis() * angular.getNlm());

    std::string filename = std::string("misc") + "/block_norms.txt";

    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(15);

    if (!outFile)
    {
        std::cerr << "Error opening file: " << filename << '\n';
    }

    if (getProjOut())
    {
        projectOutBoundStates(finalState,S,tise,angular,Basis);
    }

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        int start = blockIdx * Basis.getNbasis();

        IndexSet is{PETSC_COMM_SELF,Basis.getNbasis(), start, 1};

        Vector blockVector{};
        VecGetSubVector(finalState.get(), is.get(), &blockVector.get());

        auto normVal = innerProduct(blockVector,S,blockVector);

        VecRestoreSubVector(finalState.get(), is.get(), &blockVector.get());

        outFile << blockIdx << " " << normVal.real() << " " << normVal.imag() << '\n';
    }

    outFile.close();
}

void Observables::projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const Basis& Basis)
{
    PetscHDF5 viewer(PETSC_COMM_SELF,tise.getOutputPath(),FILE_MODE_READ);

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockIdx);


        int start = blockIdx * Basis.getNbasis();
        
        IndexSet is{PETSC_COMM_SELF,Basis.getNbasis(),start,1};

        auto blockVector = Vector{};

        VecGetSubVector(finalState.get(),is.get(),&blockVector.get());


        for (int n = 0; n <= tise.getNmax(); ++n)
        {
            std::string groupName = "eigenvectors";
            std::string vectorName = std::string("psi_") + std::to_string(n) + std::string("_") + std::to_string(l);

            std::string datasetName = groupName + "/" + vectorName;

            PetscBool hasDataset{};
            PetscViewerHDF5HasDataset(viewer.get(),datasetName.c_str(),&hasDataset);

            if (hasDataset)
            {
                Vector tiseState = viewer.loadVector(groupName,vectorName,Basis.getNbasis());

                PetscScalar prodVal = innerProduct(tiseState,S,blockVector);

                blockVector.AXPY(-prodVal, tiseState);
            }


        }

        VecRestoreSubVector(finalState.get(),is.get(),&blockVector.get());
    }
}

CoulombWave Observables::computeCoulombWave(double E, int l, const Box& box, const Atom& atom)
{
    // Unpack necessary variables
    int Nr = box.getNr();
    double dr = box.getGridSpacing();

    // Initialize empty vector to store wave
    std::vector<double> wave(Nr, 0.0);

    // Compute some relevant values
    double dr2 = dr * dr;
    double k = std::sqrt(2.0 * E);
    int lterm = l * (l + 1);

    // Bootstrap the Numerov method
    wave[0] = 0.0;
    wave[1] = 1.0;

    // Numerov method to populate wave
    for (int rIdx = 2; rIdx < Nr; ++rIdx) 
    {   
        // Compute next step
        double rVal = rIdx * dr;
        double term = dr2 * (lterm/(rVal*rVal) + 2.0*atom(rVal).real() - 2.0*E);
        wave[rIdx] = wave[rIdx - 1] * (term + 2.0) - wave[rIdx - 2];

        // Check to ensure wave hasnt blown up, and normalize if it has
        if (std::abs(wave[rIdx]) > 1E10) 
        {
            double maxMag = std::abs(*std::max_element(wave.begin(), wave.end(), 
                [](double a, double b) { return std::abs(a) < std::abs(b); }));

            for (auto& value : wave)
            {
                value /= maxMag;
            }
        }
    }

    // Compute values to normalize
    double rEnd = (Nr - 2) * dr;
    double waveEnd = wave[Nr - 2];
    double dwaveEnd = (wave[Nr - 1] - wave[Nr - 3]) / (2.0 * dr);
    
    // Normalize
    double denom = k + 1.0/(k * rEnd);
    double termPsi = std::abs(waveEnd) * std::abs(waveEnd);
    double termDer = std::abs(dwaveEnd/denom) * std::abs(dwaveEnd/denom);
    double normVal = std::sqrt(termPsi + termDer);

    if (normVal > 0.0) 
    {
        for (auto& value : wave)
        {
            value /= normVal;
        }
    }

    

    // Compute values to compute phase
    std::complex<double> numerator(0.0, waveEnd);
    numerator += dwaveEnd / denom;

    const double scale = 2.0 * k * rEnd;
    const std::complex<double> denomC = std::exp(std::complex<double>(0.0, 1.0/k) * std::log(scale));
    const std::complex<double> fraction = numerator / denomC;
    const double phase = std::arg(fraction) - k * rEnd + l * M_PI/2.0;

    return CoulombWave{wave,phase};
}

std::vector<std::complex<double>> Observables::expandState(const Vector& state, const Box& box, const Angular& angular, const Basis& basis)
{   
    int Nr = box.getNr();
    double dr = box.getGridSpacing();
    int nlm = angular.getNlm();
    int nbasis = basis.getNbasis();
    int degree = basis.getDegree();

    // Unpacked the petsc vector into a C-style array
    const std::complex<double>* stateArray;
    VecGetArrayRead(state.get(), reinterpret_cast<const PetscScalar**>(&stateArray));
    
    // Allocate vector to hold position space state
    std::vector<std::complex<double>> expanded_state(Nr * nlm);


    // Evaluate all bspline basis functions in their nonzero intervals
    for (int bsplineIdx = 0; bsplineIdx < nbasis; ++bsplineIdx)
    {   
        // Get the start and end of the interval where the bspline is nonzero
        std::complex<double> start = basis.getKnots()[bsplineIdx];
        std::complex<double> end = basis.getKnots()[bsplineIdx+degree+1];

        // Initialize vectors to store nonzero bspline values and the position indices where they occur
        std::vector<std::complex<double>> bsplineEval{};
        std::vector<int> bsplineEvalIndices{};

        // Loop over all grid points
        for (int rIdx = 0; rIdx < Nr; ++rIdx)
        {   
            // Compute position
            double r = rIdx*dr;

            // See if position is in nonzero interval for this spline, if so evaluate
            if (r >= start.real() && r < end.real())
            {
                std::complex<double> bsplineVal = BSplines::B(bsplineIdx,degree,r,basis.getKnots());
                bsplineEval.push_back(bsplineVal);
                bsplineEvalIndices.push_back(rIdx);
            }
        }

        // Loop over each block 
        for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
        {   
            // Get the bspline coefficients
            std::complex<double> coeff = stateArray[blockIdx*nbasis + bsplineIdx];

            // Loop over all grid points and add contribution to the expanded state for this block
            for (size_t rSubIdx = 0; rSubIdx < bsplineEval.size(); ++rSubIdx)
            {
                expanded_state[blockIdx*Nr + bsplineEvalIndices[rSubIdx]] += coeff*bsplineEval[rSubIdx];
            }
        }
    }
    return  expanded_state;
}

std::pair<std::map<lm_pair,std::vector<std::complex<double>>>,std::map<std::pair<double, int>,double>> Observables::computePartialSpectra(const std::vector<std::complex<double>>& expanded_state,const Angular& angular, const Atom& atom,const Box& box)
{
    std::map<lm_pair,std::vector<std::complex<double>>> partialSpectra;

    std::map<std::pair<double, int>,double> phases;

    int nlm = angular.getNlm();
    int Ne = getNe();
    double Emin = getEmin();
    int Nr = box.getNr();
    double dr = box.getGridSpacing();

    // Allocate space for each partial spectra
    for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
    {
        std::pair<int,int> lmPair = angular.getBlockMap().at(blockIdx);
        int l = lmPair.first;
        int m = lmPair.second;
        partialSpectra[std::make_pair(l, m)].reserve(Ne); 
    }


    for (int EIdx = 1; EIdx <= Ne; ++EIdx)
    {
        for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
        {
            std::pair<int,int> lm_pair = angular.getBlockMap().at(blockIdx);
            int l = lm_pair.first;
            int m = lm_pair.second;
            CoulombWave coulombResult = computeCoulombWave(EIdx*Emin, l, box,atom);
            
            phases[std::make_pair(EIdx*Emin,l)] = coulombResult.phase;


            auto start = expanded_state.begin() + Nr*blockIdx;  
            auto end = expanded_state.begin() + Nr*(blockIdx+1);    

            // Extract subvector
            std::vector<std::complex<double>> blockVector(start, end);  

            // We want to form the position space inner product aka integrate. 
            // So first we to a pointwise mult. Since coulomb wave is real no need for complex conjugation

            for (int rIdx = 0; rIdx < Nr; ++rIdx)
            {
                blockVector[rIdx] = blockVector[rIdx] * coulombResult.wave[rIdx];
            }

            std::complex<double> I = SimpsonsMethod(blockVector,dr);  
            partialSpectra[std::make_pair(l,m)].push_back(I);
        }
    }
    return std::make_pair(partialSpectra,phases);
}

void Observables::computeAngleIntegrated(const std::map<lm_pair,std::vector<std::complex<double>>>& partialSpectra,const Angular& angular)
{   
    int Ne = getNe();
    int nlm = angular.getNlm();
    double Emin = getEmin();

    std::ofstream pesFiles("misc/pes.txt", std::ios::app);
    std::vector<std::complex<double>> pes(Ne);
    
    for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
    {   
        lm_pair lm_pair = angular.getBlockMap().at(blockIdx);
        int l = lm_pair.first;
        int m = lm_pair.second;

        std::vector<std::complex<double>> partialSpectrum = partialSpectra.at(std::make_pair(l,m));

        std::vector<std::complex<double>> magSq(Ne);
        for (size_t idx = 0; idx < partialSpectrum.size(); ++idx)
        {
            magSq[idx] = partialSpectrum[idx] * std::conj(partialSpectrum[idx]);
        }

        for (size_t idx = 0; idx < magSq.size(); ++idx)
        {
            pes[idx] += magSq[idx];
        }
    }
    

    for (size_t idx = 1; idx < pes.size(); ++idx)
    {   
        std::complex<double> val = pes[idx];
        val /= ((2*M_PI)*(2*M_PI)*(2*M_PI));
        pesFiles << idx*Emin << " " << val.real() << " " << "\n";
    }

    pesFiles.close();
}

void Observables::computeAngleResolved(const std::map<lm_pair,std::vector<std::complex<double>>>& partialSpectra,std::map<std::pair<double, int>,double> phases)
{   
    int Ne = getNe();
    double Emin = getEmin();

        std::ofstream padFiles("misc/pad.txt", std::ios::app);
        std::vector<double> theta_range;
        std::vector<double> phi_range;

        if (getSlice() == "XZ")
        {
            for (double theta = 0; theta <= M_PI; theta += 0.01) 
            {
                theta_range.push_back(theta);
            }

            phi_range  = {0.0,M_PI};
        }
        if (getSlice() == "XY")
        {
            theta_range = {M_PI/ 2.0};

            for (double phi = 0; phi < 2.0*M_PI; phi += 0.01) 
            {
                phi_range.push_back(phi);
            }

        }
        for (int EIdx = 1; EIdx <= Ne; ++EIdx)
        {   

            double E = EIdx*Emin;
            double k = std::sqrt(2.0*E);

            for (auto& theta : theta_range)
            {
                for (auto& phi : phi_range)
                {   
                    std::complex<double> pad_amplitude {};
                    for (const auto& pair: partialSpectra)
                    {
                        int l = pair.first.first;
                        int m = pair.first.second;

                        std::complex<double> sph_term = compute_Ylm(l,m,theta,phi);

                        double partial_phase = phases[std::make_pair(E,l)];
                        std::complex<double> partial_amplitude = pair.second[EIdx - 1];

                        std::complex<double> phase_factor = std::exp(std::complex<double>(0.0,3*l*M_PI/2.0 + partial_phase));

                        pad_amplitude += sph_term*phase_factor*partial_amplitude;

                    }

                    double pad_prob = std::norm(pad_amplitude);
                    pad_prob/=((2*M_PI)*(2*M_PI)*(2*M_PI));
                    pad_prob/=k;

                    padFiles << E << " " << theta << " " << phi << " " << pad_prob << "\n";
                }
            }
        }
}

void Observables::computePhotoelectronSpectrum(int rank,const TISE& tise, const TDSE& tdse, const Angular& angular, const Basis& basis, const Box& box, const Atom& atom)
{
    if (rank != 0)
    {
        return;
    }

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
    RadialMatrix::populateRadialMatrix(RadialMatrixType::S,S,basis,false);

    PetscHDF5 viewer{PETSC_COMM_SELF,tdse.getOutputPath(),FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,basis.getNbasis() * angular.getNlm());

    projectOutBoundStates(finalState,S,tise,angular,basis);

    std::vector<std::complex<double>> expandedState = expandState(finalState,box,angular,basis);


    std::map<lm_pair,std::vector<std::complex<double>>> partialSpectra;
    std::map<std::pair<double, int>,double> phases;

    std::tie(partialSpectra,phases) = computePartialSpectra(expandedState,angular,atom,box);

    computeAngleIntegrated(partialSpectra,angular);
    computeAngleResolved(partialSpectra,phases);
}

void Observables::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        std::cout << "Observables Configuration: " << "\n\n";
        std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
        std::cout << "Block :  projOutBound: " << getProjOut() << "\n\n";
    }
}