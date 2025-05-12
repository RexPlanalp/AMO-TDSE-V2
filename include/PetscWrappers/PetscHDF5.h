#include "petscviewerhdf5.h"
#include "PetscWrappers/PetscVec.h"

class PetscHDF5
{
    public:
        

        PetscHDF5(MPI_Comm comm, const std::string& fileName, PetscFileMode fileMode) 
        {
            PetscViewerHDF5Open(comm,fileName.c_str(),fileMode,&viewer);
        }

        void saveVector(const std::string& groupName, const std::string& vectorName, Vector& vector)
        {
            PetscObjectSetName(PetscObject(vector.get()), vectorName.c_str());

            PetscViewerHDF5PushGroup(viewer,groupName.c_str());
            VecView(vector.get(),viewer);
            PetscViewerHDF5PopGroup(viewer);
        }

        void saveValue(const std::string& groupName, const std::string& valueName, PetscScalar value)
        {
            auto temp = Vector{PETSC_COMM_WORLD, PETSC_DETERMINE, 1};
            temp.setValue(0,value);
            temp.assemble();

            saveVector(groupName,valueName,temp);
        }

        Vector loadVector(const std::string& groupName, const std::string& vectorName, PetscInt size)
        {
            auto temp = Vector{PETSC_COMM_SELF,PETSC_DECIDE,size};
            auto totalName = groupName + "/" + vectorName;
            
            PetscObjectSetName(PetscObject(temp.get()), totalName.c_str());

            VecLoad(temp.get(), viewer);

            return temp;
        }

    private:
        PetscViewer viewer{nullptr};
        std::string fileName{};
   


};