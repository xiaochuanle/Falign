#include "../../corelib/hbn_aux.h"
#include "command_list.hpp"

#include <cstring>

using namespace std;

int main(int argc, char* argv[])
{
    Dip3dCommandList cmds;
    int r = cmds.run_cmd(argc, argv);
    if (r) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s COMMANDS\n", argv[0]);
        cmds.dump_cmds();
        return 1;
    }

    return 0;
}