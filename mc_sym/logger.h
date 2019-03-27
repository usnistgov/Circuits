#pragma once

#include <iostream>


#define LOG_LEVEL_SILENT	0
#define LOG_LEVEL_ERROR		1
#define LOG_LEVEL_WARN		2
#define LOG_LEVEL_INFO		3
#define LOG_LEVEL_DEBUG		4
#define LOG_LEVEL_TRACE		5
#define LOG_LEVEL_VERBOSE	0xff


#ifndef LOG_LEVEL
	#define LOG_LEVEL	LOG_LEVEL_INFO
#endif


struct empty_logger{};

struct console_logger {};

struct logger_endl {};

template<typename T>
empty_logger& operator<<(empty_logger &v, const T& t) {

	return v;
}

empty_logger& operator<<(empty_logger &v, const char* t);

empty_logger& operator<<(empty_logger &v, const logger_endl&);

template<typename T>
console_logger& operator<<(console_logger &v, const T& t) {
	
	std::cout << t;
	return v;
}

console_logger& operator<<(console_logger &v, const char* t);

console_logger& operator<<(console_logger &v, const logger_endl&);

extern empty_logger		g_empty_logger;
extern console_logger	g_console_logger;
extern logger_endl		g_logger_endl;

#define LOG_CONSOLE		g_console_logger
#define LOG_NOTHING		g_empty_logger
#define LOG_ENDL		g_logger_endl
#define LOG_CONSOLE_OMP	#pragma omp critical LOG_CONSOLE

#if LOG_LEVEL_ERROR <= LOG_LEVEL
	#define LOG_ERROR			LOG_CONSOLE
	#define LOG_ERROR_OMP		LOG_CONSOLE_OMP
	#define LOG_ERROR_IF(cond)	if(cond) LOG_ERROR
#else
	#define LOG_ERROR			LOG_NOTHING 
	#define LOG_ERROR_OMP		LOG_NOTHING 
	#define LOG_ERROR_IF(cond)	LOG_NOTHING 
#endif

#if LOG_LEVEL_WARN <= LOG_LEVEL
	#define LOG_WARN			LOG_CONSOLE
	#define LOG_WARN_OMP		LOG_CONSOLE_OMP
	#define LOG_WARN_IF(cond)	if(cond) LOG_WARN
#else
	#define LOG_WARN			LOG_NOTHING 
	#define LOG_WARN_OMP		LOG_NOTHING 
	#define LOG_WARN_IF(cond)	LOG_NOTHING 
#endif

#if LOG_LEVEL_INFO <= LOG_LEVEL
	#define LOG_INFO			LOG_CONSOLE
	#define LOG_INFO_OMP		LOG_CONSOLE_OMP
	#define LOG_INFO_IF(cond)	if(cond) LOG_INFO
#else
	#define LOG_INFO			LOG_NOTHING 
	#define LOG_INFO_OMP		LOG_NOTHING 
	#define LOG_INFO_IF(cond)	LOG_NOTHING 
#endif

#if LOG_LEVEL_DEBUG <= LOG_LEVEL
	#define LOG_DEBUG			LOG_CONSOLE
	#define LOG_DEBUG_OMP		LOG_CONSOLE_OMP
	#define LOG_DEBUG_IF(cond)	if(cond) LOG_DEBUG
#else
	#define LOG_DEBUG			LOG_NOTHING 
	#define LOG_DEBUG_OMP		LOG_NOTHING 
	#define LOG_DEBUG_IF(cond)	LOG_NOTHING 
#endif

#if LOG_LEVEL_TRACE <= LOG_LEVEL
	#define LOG_TRACE			LOG_CONSOLE
	#define LOG_TRACE_OMP		LOG_CONSOLE_OMP
	#define LOG_TRACE_IF(cond)	if(cond) LOG_TRACE
#else
	#define LOG_TRACE			LOG_NOTHING
	#define LOG_TRACE_OMP		LOG_NOTHING
	#define LOG_TRACE_IF(cond)	LOG_NOTHING
#endif



//
// Assertion and Verification 
//

#define LOG_VERIFY(expr) if(!(expr)) { LOG_CONSOLE << "assertion '" << #expr << "' failed" << '\n';}

#define LOG_VERIFY_IF(cond, expr) if((cond) && !(expr)) { LOG_CONSOLE << "assertion '" << #expr << "' failed" << '\n';}

#define LOG_VERIFY_VERBOSE(expr)  if(expr) {	LOG_CONSOLE << "assertion '" << #expr << "' verified" << LOG_ENDL; } else { \
												LOG_CONSOLE << "assertion '" << #expr << "' failed"   << LOG_ENDL; }

#define LOG_VERIFY_EQ(expr1, expr2) { auto var1 = expr1; auto var2 = expr2; \
										  if(var1 != var2) { \
											LOG_CONSOLE << "assertion '" << #expr1 << " == " << #expr2 << "' failed" << LOG_ENDL;\
											LOG_CONSOLE << #expr1 << " = " << var1 << LOG_ENDL;\
											LOG_CONSOLE << #expr2 << " = " << var2 << '\n' << LOG_ENDL; } \
										}


#ifdef NDEBUG
	#define LOG_ASSERT(expr)			((void)0)
	#define LOG_ASSERT_IF(cond, expr)	((void)0)
	#define LOG_ASSERT_VERBOSE(expr)	((void)0)
	#define LOG_ASSERT_EQ(expr1, expr2) ((void)0)
#else
	#define LOG_ASSERT(expr)			LOG_VERIFY(expr)
	#define LOG_ASSERT_IF(cond, expr)	if(cond) { LOG_ASSERT(expr); }
	#define LOG_ASSERT_VERBOSE(expr)	LOG_VERIFY_VERBOSE(expr)
	#define LOG_ASSERT_EQ(expr1, expr2) LOG_VERIFY_EQ(expr1, expr2)
#endif
