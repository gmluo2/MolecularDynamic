// MDInforDlg.cpp : implementation file
//

#include "stdafx.h"
#include "NEIMOMD_DLL.h"
#include "MDInforDlg.h"


// CMDInforDlg dialog

IMPLEMENT_DYNAMIC(CMDInforDlg, CDialog)

CMDInforDlg::CMDInforDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMDInforDlg::IDD, pParent)
	, md_loop(0)
	, md_nsaved(0)
	, md_torun(0)
	, md_Ek(0)
	, md_Ep(0)
	, md_Etotal(0)
	, md_infor(_T(""))
{
	ncur = 0;
	int i = 0;
	for (i = 0; i < INFOR_LINES; i++) strcpy(msg[i], "\0");
}

CMDInforDlg::~CMDInforDlg()
{
}

void CMDInforDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_LOOP, md_loop);
	DDX_Text(pDX, IDC_EDIT_SAVED, md_nsaved);
	DDX_Text(pDX, IDC_EDIT_TORUN, md_torun);
	DDX_Text(pDX, IDC_EDIT_EK, md_Ek);
	DDX_Text(pDX, IDC_EDIT_POTENTIAL, md_Ep);
	DDX_Text(pDX, IDC_EDIT_ETOTAL, md_Etotal);
	DDX_Text(pDX, IDC_EDIT_INFOR, md_infor);
	DDX_Control(pDX, IDC_EDIT_INFOR, cInfor);
}


BEGIN_MESSAGE_MAP(CMDInforDlg, CDialog)
	ON_WM_CREATE()
	ON_MESSAGE(WM_DIALOG_SHOW, OnUpdateDialog)
	ON_MESSAGE(WM_DIALOG_SHOW_LOOP, OnShowLoop)
	ON_MESSAGE(WM_DIALOG_SHOW_INFOR, OnShowInfor)
	ON_WM_DESTROY()
END_MESSAGE_MAP()


// CMDInforDlg message handlers

void CMDInforDlg::append_msg(char *msg) {
	strcpy(this->msg[ncur], msg);
	strcat(this->msg[ncur], "\r\n");
	ncur++;
	if (ncur >= INFOR_LINES) ncur = 0;
	int i = 0, n = 0;
	for (i = 0; i < INFOR_LINES; i++) {
		n = i + ncur;
		while (n >= INFOR_LINES) n -= INFOR_LINES;
		if (i == 0) md_infor = CString(this->msg[n]);
		else md_infor += CString(this->msg[n]);
	}
}

int CMDInforDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  Add your specialized creation code here

	return 0;
}

LRESULT CMDInforDlg::OnUpdateDialog(WPARAM wParam, LPARAM lParam) {
	UpdateData(FALSE);
	return 1;
}

void CMDInforDlg::OnDestroy()
{
	CDialog::OnDestroy();

	// TODO: Add your message handler code here
}

LRESULT CMDInforDlg::OnShowLoop(WPARAM wParam, LPARAM lParam) {
	CDataExchange dx(this, FALSE), *pDX = &dx;
	//CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_LOOP, md_loop);
	DDX_Text(pDX, IDC_EDIT_SAVED, md_nsaved);
	DDX_Text(pDX, IDC_EDIT_TORUN, md_torun);
	DDX_Text(pDX, IDC_EDIT_EK, md_Ek);
	DDX_Text(pDX, IDC_EDIT_POTENTIAL, md_Ep);
	DDX_Text(pDX, IDC_EDIT_ETOTAL, md_Etotal);

	return 1;
}

LRESULT CMDInforDlg::OnShowInfor(WPARAM wParam, LPARAM lParam) {
	CDataExchange dx(this, FALSE), *pDX = &dx;
	//CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_INFOR, md_infor);
	cInfor.LineScroll(2);
	return 1;
}
