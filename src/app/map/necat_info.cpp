#include "necat_info.hpp"

#include <sys/utsname.h>
#include <string>
#include <thread>

#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

using namespace std;

extern "C"
size_t getMemorySizeBytes();

void dump_necat_info(FILE* out)
{
    struct utsname _os_info_buf;
    struct utsname* os_info = nullptr;
    if (uname(&_os_info_buf) == 0) os_info = &_os_info_buf;

    size_t sys_mem_bytes = getMemorySizeBytes();
    string sys_mem = NStr::UInt8ToString_DataSize(sys_mem_bytes);

    int cpu_threads = thread::hardware_concurrency();

    fprintf(out, "\n");
    fprintf(out, "PROGRAM:\n");
    fprintf(out, "  Name:           %s\n", HBN_PACKAGE_NAME);
    fprintf(out, "  Version:        %s\n", HBN_PACKAGE_VERSION);
    fprintf(out, "  Description:    Alignment toolkit for long noisy chromosome conformation capture (3C) reads\n");
    fprintf(out, "  Contact:        chenying2016@gmail.com\n");

    fprintf(out, "\n");
    fprintf(out, "SYSTEM:\n");
    if (os_info) {
    fprintf(out, "  Computer:       %s\n", os_info->nodename);
    fprintf(out, "  Name:           %s\n", os_info->sysname);
    fprintf(out, "  Release:        %s\n", os_info->release);
    fprintf(out, "  Version:        %s\n", os_info->version);
    fprintf(out, "  Machine:        %s\n", os_info->machine);
    }
    fprintf(out, "  CPU threads:    %d\n", cpu_threads);
    fprintf(out, "  RAM:            %s\n", sys_mem.c_str());
    fprintf(out, "\n");
}

extern "C" size_t getPeakRSS();

HbnRunningInfo::~HbnRunningInfo() {
    gettimeofday(&M_end, NULL);
    double dur = hbn_time_diff(&M_begin, &M_end);
    size_t peak_ram = getPeakRSS();
    string size = NStr::UInt8ToString_DataSize(peak_ram);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "%s Wallclock time: %.2f seconds.\n", HBN_PACKAGE_NAME, dur);
    fprintf(stderr, "%s Peak RAM usage: %s\n", HBN_PACKAGE_NAME, size.c_str());
}