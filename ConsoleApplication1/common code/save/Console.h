#ifndef __CONSOLE_H__
#define __CONSOLE_H__


#define PLATFORM_ENABLE_CONSOLE

#endif //__CONSOLE_H__

/* How to use it :

 1. #include "Console.h"   in the app header file.
 2. Add
			#ifdef PLATFORM_ENABLE_CONSOLE
			extern void RedirectIOToConsole(const wchar_t title[]);
			#endif // PLATFORM_ENABLE_CONSOLE
    
 3. Add
			#ifdef 	PLATFORM_ENABLE_CONSOLE
			RedirectIOToConsole(_T("Console"));
			#endif //PLATFORM_ENABLE_CONSOLE
*/