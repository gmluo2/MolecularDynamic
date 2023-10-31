#pragma once
#include "afxwin.h"


// CMDInforDlg dialog
#define INFOR_LINES  2

class CMDInforDlg : public CDialog
{
	DECLARE_DYNAMIC(CMDInforDlg)

public:
	CMDInforDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CMDInforDlg();

// Dialog Data
	enum { IDD = IDD_DIALOGBAR_MDINFOR };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

private:
	char msg[INFOR_LINES][2560];
	int ncur;
public:
	void append_msg(char *msg);
public:
	// current loop ... 
	int md_loop;
public:
	int md_nsaved;
public:
	int md_torun;
public:
	float md_Ek;
public:
	float md_Ep;
public:
	float md_Etotal;
public:
	CString md_infor;
public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
public:
	afx_msg void OnDestroy();
public:
	CEdit cInfor;
public:
	afx_msg LRESULT OnUpdateDialog(WPARAM wParam, LPARAM lParam);
public:
	afx_msg LRESULT OnShowLoop(WPARAM wParam, LPARAM lParam);
public:
	afx_msg LRESULT OnShowInfor(WPARAM wParam, LPARAM lParam);
};
