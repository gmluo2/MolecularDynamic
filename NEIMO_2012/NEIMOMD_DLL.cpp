// NEIMOMD_DLL.cpp : Defines the initialization routines for the DLL.
//

#include "stdafx.h"
#include "TAMD\def.h"
#include "NEIMOMD_DLL.h"
#include "MDInforDlg.h"
#include "MolOperateDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

extern HANDLE hWaitGetMolStructEvent;
extern HANDLE hWaitEvent;
CMDInforDlg *mdInforDlg = NULL;
CMolOperateDlg *mMolOperateDlg = NULL;
extern CWnd *pMainFrame;

//
//TODO: If this DLL is dynamically linked against the MFC DLLs,
//		any functions exported from this DLL which call into
//		MFC must have the AFX_MANAGE_STATE macro added at the
//		very beginning of the function.
//
//		For example:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// normal function body here
//		}
//
//		It is very important that this macro appear in each
//		function, prior to any calls into MFC.  This means that
//		it must appear as the first statement within the 
//		function, even before any object variable declarations
//		as their constructors may generate calls into the MFC
//		DLL.
//
//		Please see MFC Technical Notes 33 and 58 for additional
//		details.
//


// CNEIMOMD_DLLApp

BEGIN_MESSAGE_MAP(CNEIMOMD_DLLApp, CWinApp)
END_MESSAGE_MAP()


// CNEIMOMD_DLLApp construction

CNEIMOMD_DLLApp::CNEIMOMD_DLLApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only CNEIMOMD_DLLApp object

CNEIMOMD_DLLApp theApp;


// CNEIMOMD_DLLApp initialization

extern "C" bool init_env();
extern "C" void clear_env();
extern char errmsg[2560];
extern void show_msg(char *msg, bool endline);

BOOL CNEIMOMD_DLLApp::InitInstance()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	CWinApp::InitInstance();
	hWaitGetMolStructEvent = CreateEvent(NULL, FALSE, TRUE, NULL);
	hWaitEvent = CreateEvent(NULL, FALSE, FALSE, NULL);

	mMolOperateDlg = new CMolOperateDlg;
	mMolOperateDlg->Create(IDD_DIALOG_CONSTRUCT_MMOLECULE, NULL);
	mMolOperateDlg->m_MolType1.SetCurSel(0);
	mMolOperateDlg->m_MolType2.SetCurSel(0);
	mMolOperateDlg->show_dialog();

	if (sizeof(double) != SIZE_DOUBLE || sizeof(float) != SIZE_FLOAT ||
		sizeof(int) != SIZE_INT || sizeof(char) != SIZE_CHAR) {
		show_msg("Warning: incorrect definition of SIZE_DOUBLE, and related definitions", false);
		{
			char msg[256] = "\0";
			sprintf(msg, "variable size: char -- %d; int -- %d; float -- %d; double -- %d", sizeof(char), sizeof(int), sizeof(float), sizeof(double));
			show_msg(msg, false);
		}
	}

	if (!init_env()) {
		//show_msg(::errmsg, true);
		//return FALSE;
	}
	/*
	extern void test_GaussRand(float width, int max);
	int nmax = 100000;
	test_GaussRand(1, nmax);
	*/
	return TRUE;
}

int CNEIMOMD_DLLApp::ExitInstance()
{
	// TODO: Add your specialized code here and/or call the base class
	if (mdInforDlg != NULL) {
		if (mdInforDlg->m_hWnd != NULL) SendMessage(mdInforDlg->m_hWnd, WM_DESTROY, 0, 0);
		delete mdInforDlg; mdInforDlg = NULL;
	}
	if (mMolOperateDlg != NULL) {
		if (mMolOperateDlg->m_hWnd != NULL) SendMessage(mMolOperateDlg->m_hWnd, WM_DESTROY, 0, 0);
		delete mMolOperateDlg; mMolOperateDlg = NULL;
	}

	CloseHandle(hWaitGetMolStructEvent);
	CloseHandle(hWaitEvent);

	clear_env();

	return CWinApp::ExitInstance();
}

void create_infor_dialog() {
	if (mdInforDlg != NULL) delete mdInforDlg;
	mdInforDlg = new CMDInforDlg;
	mdInforDlg->Create(IDD_DIALOGBAR_MDINFOR, NULL); //pMainFrame);
	mdInforDlg->ShowWindow(SW_SHOWDEFAULT); //(SW_SHOW);
}

void show_mol_operate_dialog() {
	if (!mMolOperateDlg->show) mMolOperateDlg->ShowWindow(SW_SHOW);
}
