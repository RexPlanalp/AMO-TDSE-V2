#include "misc.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"

int main()
{   

    nlohmann::json J = loadJson("input.json");
    std::array<int,3> components{1,0,1};
    std::array<int,3> initial_state{1,1,0};
    int rank = 0;

    Box box{J};
    Angular angular{J};
    Laser laser{J};

    angular.buildMaps(components,initial_state);
    angular.dumpTo("misc",rank);
    laser.dumpTo("misc",rank);

  

    return 0;
}