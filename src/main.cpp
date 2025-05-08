#include "misc.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"
#include "BSpline.h"
#include "Atom.h"

int main()
{   

    nlohmann::json J = loadJson("input.json");
    std::array<int,3> initial_state{1,1,0};
    int rank = 0;

    Box box{J};
    Angular angular{J};
    Laser laser{J};
    Atom atom{J};
    BSpline bspline{J,box};


    angular.buildMaps(laser,initial_state);



    bspline.dumpTo(box,"misc",rank);
    angular.dumpTo("misc",rank);
    laser.dumpTo("misc",rank);

  

    return 0;
}