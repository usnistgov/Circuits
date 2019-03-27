#include "platform.h"

static PlatformInfo g_this_platform;

const PlatformInfo &get_platform_info() {
	return g_this_platform;
}

void PlatformInfo::populate()
{


#ifdef _WIN32

	os = OperatingSystem::Windows;

	#ifdef _WIN64
		arch = Architecture::x64;
	#else
		arch = Architecture::x86;
	#endif

#elif defined( __APPLE__) && defined(__MACH__)

	os = OperatingSystem::Mac;

#elif defined(__linux__)

	os = OperatingSystem::Linux;

#elif defined(__unix__)

	os = OperatingSystem::Unix;

#endif


#if defined(_MSC_VER)
	compiler = Compiler::MSC;
#elif defined(__INTEL_COMPILER)
	compiler = Compiler::Intel;
#elif defined(__GNUC__)
	compiler = Compiler::GCC;
#elif defined(__clang__)
	compiler = Compiler::Clang;
#endif
}
