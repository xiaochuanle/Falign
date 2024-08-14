#ifndef __NECAT_SYS_INFO_HPP
#define __NECAT_SYS_INFO_HPP

#include <cstdio>

#include "../../corelib/hbn_aux.h"

#define HBN_PACKAGE                         1
#define HBN_PACKAGE_NAME                    "falign"
#define HBN_PACKAGE_VERSION_MAJOR           2
#define HBN_PACKAGE_VERSION_MINOR           0
#define HBN_PACKAGE_VERSION_PATCH           0
#define HBN_PACKAGE_CONFIG                  ""

#define HBN_PACKAGE_VERSION_STRINGIFY(x)    #x
#define HBN_PACKAGE_VERSION_COMPOSE_STR(a, b, c)    \
    HBN_PACKAGE_VERSION_STRINGIFY(a) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(b) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(c)

#define HBN_PACKAGE_VERSION             \
    HBN_PACKAGE_VERSION_COMPOSE_STR     \
    (                                   \
        HBN_PACKAGE_VERSION_MAJOR,      \
        HBN_PACKAGE_VERSION_MINOR,      \
        HBN_PACKAGE_VERSION_PATCH       \
    )

void dump_necat_info(FILE* out = stderr);

struct HbnRunningInfo
{
public:
    HbnRunningInfo() {
        dump_necat_info();
        gettimeofday(&M_begin, NULL);
    }

    ~HbnRunningInfo();
    
private:
    struct timeval M_begin;
    struct timeval M_end;
};

#endif // __NECAT_SYS_INFO_HPP