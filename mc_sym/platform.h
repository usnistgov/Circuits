#pragma once


//
// TODO: configure instructions supported by a CPU
#define ENABLE_HARDWARE_POPCNT

// _BitScanReverse


enum class Architecture {Unknown, x86, x64};

enum class OperatingSystem {Unknown, Linux, Mac, Unix, Windows };

enum class Compiler {Unknown, GCC, Clang, Intel, MSC};

struct PlatformInfo {

	PlatformInfo() {
		populate();
	}

	Architecture	arch		= Architecture::Unknown;
	OperatingSystem os			= OperatingSystem::Unknown;
	Compiler		compiler	= Compiler::Unknown;

private:
	void populate();
};

const PlatformInfo& get_platform_info();
