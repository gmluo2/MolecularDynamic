// NEIMOMD_DLL.h : main header file for the NEIMOMD_DLL DLL
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CNEIMOMD_DLLApp
// See NEIMOMD_DLL.cpp for the implementation of this class
//

class CNEIMOMD_DLLApp : public CWinApp
{
public:
	CNEIMOMD_DLLApp();

// Overrides
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
public:
	virtual int ExitInstance();
};
