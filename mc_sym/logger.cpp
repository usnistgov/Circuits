#include "logger.h"

empty_logger	g_empty_logger;
console_logger	g_console_logger;
logger_endl		g_logger_endl;

empty_logger& operator<<(empty_logger &v, const logger_endl&) {

	return v;
}

empty_logger& operator<<(empty_logger &v, const char* t) {

	return v;
}

console_logger& operator<<(console_logger &v, const logger_endl&) {

	std::cout << std::endl;
	return v;
}

console_logger& operator<<(console_logger &v, const char* t) {

	std::cout << t;
	return v;
}
