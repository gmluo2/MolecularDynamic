#pragma once
#include "afxwin.h"


// CMolOperateDlg dialog

class CMolOperateDlg : public CDialog
{
	DECLARE_DYNAMIC(CMolOperateDlg)

public:
	bool show;
	CMolOperateDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CMolOperateDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_CONSTRUCT_MMOLECULE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CButton mCheck_MolName;
public:
	CEdit mEdit_MolName;
public:
	CButton mButton_MolCreate;
public:
	CButton mCheck_ReadMol;
public:
	CEdit mEdit_MolFilename;
public:
	CButton mButton_FileBrowse;
public:
	CButton mButton_ReadMol;
public:
	CComboBox mCommand_Box;
public:
	BOOL bMolName;
public:
	BOOL bReadMolFromFile;
public:
	CString mMolFilename;
public:
	CString mMolName;
public:
	CString mCommand;
public:
	afx_msg void OnBnClickedCheckMolname();
public:
	afx_msg void OnBnClickedCheckReadfromfile();
public:
	afx_msg void OnBnClickedButton2();
public:
	afx_msg void OnBnClickedButtonCreateMolname();
public:
	afx_msg void OnBnClickedButtonReadmolFromfile();
public:
	CButton mCheck_RunCommand;
public:
	BOOL bRunCommand;
public:
	CButton mButton_RunCommand;
public:
	afx_msg void OnBnClickedButton3();
public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
public:
	void show_mol_name();
	void show_read_mol();
	void show_run_command();
	void show_multi_mol();
	void show_dialog();
public:
	afx_msg void OnBnClickedCheckRuncommand();
public:
	void check_mol_name();
	void check_read_from_file();
	void check_read_multi_mol();
	void check_run_command();
public:
	CEdit mEdit_Base;
public:
	int nBase;
public:
	CEdit mEdit_TAngle;
public:
	CString mTAngle;
public:
	afx_msg void OnClose();
public:
	afx_msg void OnBnClickedButtonClearMonitor();
public:
	CButton mCheck_MultiMol;
public:
	BOOL bMultiMol;
public:
	CEdit mEdit_MultiMolFname;
public:
	CString mMultiMolFname;
public:
	CButton mButtom_BrowseMultiMolFname;
public:
	CButton mButton_ReadMultiMol;
public:
	afx_msg void OnBnClickedCheckReadmultiMol();
public:
	afx_msg void OnBnClickedButtonBrowseMultiMolFname();
public:
	afx_msg void OnBnClickedButtonReadMultiMolFromfile();
public:
	CComboBox m_MolType1;
public:
	CComboBox m_MolType2;
};
