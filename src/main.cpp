#include "misc.h"
#include "Box.h"
#include "Angular.h"

int main()
{   

    nlohmann::json J = loadJson("input.json");

    Box box{J};
    Angular angular{J};

    std::array<int,3> components{1,0,1};
    std::array<int,3> initial_state{1,1,0};

    angular.buildMaps(components,initial_state);
    angular.dumpTo("misc");


  

    return 0;
}