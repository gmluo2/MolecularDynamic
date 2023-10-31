// MolOperateDlg.cpp : implementation file
//
#define _CRT_SECURE_NO_DEPRECATE 1
#include "stdafx.h"
#include "NEIMOMD_DLL.h"
#include "MolOperateDlg.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "TAMD\def.h"
#include "TAMD\show.h"
#include "TAMD\ranlib.h"
#include "TAMD\vector.h"
#include "TAMD\bound.h"
#include "TAMD\Matrix.h"
#include "TAMD\nhc.h"
#include "TAMD\EwaldSumVar.h"
#include "TAMD\MM.h"
#include "TAMD\cg-mm.h"
#include "TAMD\Mol.h"
#include "TAMD\CMM.h"
#include "TAMD\MD.h"
#include "TAMD\var.h"
#include "TAMD\cluster.h"
#include "TAMD\cg-cluster.h"
#include "TAMD\command.h"
#include "TAMD\Interaction.h"
#include "TAMD\Interact1.h"
#include "TAMD\md_cell.h"
#include "TAMD\mol-cell.h"

using namespace _coarse_grain_;

extern LIST<MMOLECULE> *mlist;

extern LIST<CG_MMOLECULE> *cgmlist;

extern COMMAND* pCmd;
extern COMMAND* pCurrentCmd;
extern char cmd_msg[256];
extern int nCurMol;

extern void Upper(char* str);
extern bool NextSection(char** original, char ch);
extern bool save_xyz_struct(MMOLECULE &mm, char *fname);
extern bool search(ifstream &in, char *title, char *buffer);
extern bool ConstructMoleculeStruct(MMOLECULE *m, char *fname);
extern bool ConstructCGMMStruct(CG_MMOLECULE *m, char *fname);
extern CWnd *pMainFrame;

void string2char(CString &str, char *msg) {
	int i = 0;
	for (i = 0; i < str.GetLength(); i++) msg[i] = str[i];
	msg[i] = '\0';
}

// CMolOperateDlg dialog

IMPLEMENT_DYNAMIC(CMolOperateDlg, CDialog)

CMolOperateDlg::CMolOperateDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMolOperateDlg::IDD, pParent)
	, bMolName(FALSE)
	, bReadMolFromFile(FALSE)
	, mMolFilename(_T(""))
	, mMolName(_T(""))
	, mCommand(_T(""))
	, bRunCommand(FALSE)
	, nBase(0)
	, mTAngle(_T(""))
	, bMultiMol(FALSE)
	, mMultiMolFname(_T(""))
{
	show = false;
}

CMolOperateDlg::~CMolOperateDlg()
{
}

void CMolOperateDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_CHECK_MOLNAME, mCheck_MolName);
	DDX_Control(pDX, IDC_EDIT_MOLECULE_NAME, mEdit_MolName);
	DDX_Control(pDX, IDC_BUTTON_CREATE_MOLNAME, mButton_MolCreate);
	DDX_Control(pDX, IDC_CHECK_READFROMFILE, mCheck_ReadMol);
	DDX_Control(pDX, IDC_EDIT_MOL_FNAME, mEdit_MolFilename);
	DDX_Control(pDX, IDC_BUTTON2, mButton_FileBrowse);
	DDX_Control(pDX, IDC_BUTTON_READMOL_FROMFILE, mButton_ReadMol);
	DDX_Control(pDX, IDC_COMBO_COMMAND, mCommand_Box);
	DDX_Check(pDX, IDC_CHECK_MOLNAME, bMolName);
	DDX_Check(pDX, IDC_CHECK_READFROMFILE, bReadMolFromFile);
	DDX_Check(pDX, IDC_CHECK_RUNCOMMAND, bRunCommand);
	DDX_Text(pDX, IDC_EDIT_MOL_FNAME, mMolFilename);
	DDX_Text(pDX, IDC_EDIT_MOLECULE_NAME, mMolName);
	DDX_CBString(pDX, IDC_COMBO_COMMAND, mCommand);
	DDX_Control(pDX, IDC_CHECK_RUNCOMMAND, mCheck_RunCommand);
	DDX_Control(pDX, IDC_BUTTON3, mButton_RunCommand);
	DDX_Control(pDX, IDC_EDIT_NBASE, mEdit_Base);
	DDX_Text(pDX, IDC_EDIT_NBASE, nBase);
	DDX_Control(pDX, IDC_EDIT_MOLECULE_TANGLE, mEdit_TAngle);
	DDX_Text(pDX, IDC_EDIT_MOLECULE_TANGLE, mTAngle);
	DDX_Control(pDX, IDC_CHECK_READMULTI_MOL, mCheck_MultiMol);
	DDX_Check(pDX, IDC_CHECK_READMULTI_MOL, bMultiMol);
	DDX_Control(pDX, IDC_EDIT_MULTI_MOL_FNAME, mEdit_MultiMolFname);
	DDX_Text(pDX, IDC_EDIT_MULTI_MOL_FNAME, mMultiMolFname);
	DDX_Control(pDX, IDC_BUTTON_BROWSE_MULTI_MOL_FNAME, mButtom_BrowseMultiMolFname);
	DDX_Control(pDX, IDC_BUTTON_READ_MULTI_MOL_FROMFILE, mButton_ReadMultiMol);
	DDX_Control(pDX, IDC_COMBO_MOL_TYPE1, m_MolType1);
	DDX_Control(pDX, IDC_COMBO_MOL_TYPE2, m_MolType2);
}


BEGIN_MESSAGE_MAP(CMolOperateDlg, CDialog)
	ON_BN_CLICKED(IDC_CHECK_MOLNAME, &CMolOperateDlg::OnBnClickedCheckMolname)
	ON_BN_CLICKED(IDC_CHECK_READFROMFILE, &CMolOperateDlg::OnBnClickedCheckReadfromfile)
	ON_BN_CLICKED(IDC_BUTTON2, &CMolOperateDlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON_CREATE_MOLNAME, &CMolOperateDlg::OnBnClickedButtonCreateMolname)
	ON_BN_CLICKED(IDC_BUTTON_READMOL_FROMFILE, &CMolOperateDlg::OnBnClickedButtonReadmolFromfile)
	ON_BN_CLICKED(IDC_BUTTON3, &CMolOperateDlg::OnBnClickedButton3)
	ON_WM_CREATE()
	ON_BN_CLICKED(IDC_CHECK_RUNCOMMAND, &CMolOperateDlg::OnBnClickedCheckRuncommand)
	ON_WM_CLOSE()
	ON_BN_CLICKED(IDC_BUTTON_CLEAR_MONITOR, &CMolOperateDlg::OnBnClickedButtonClearMonitor)
	ON_BN_CLICKED(IDC_CHECK_READMULTI_MOL, &CMolOperateDlg::OnBnClickedCheckReadmultiMol)
	ON_BN_CLICKED(IDC_BUTTON_BROWSE_MULTI_MOL_FNAME, &CMolOperateDlg::OnBnClickedButtonBrowseMultiMolFname)
	ON_BN_CLICKED(IDC_BUTTON_READ_MULTI_MOL_FROMFILE, &CMolOperateDlg::OnBnClickedButtonReadMultiMolFromfile)
END_MESSAGE_MAP()


void CMolOperateDlg::show_mol_name() {
	this->mEdit_MolName.EnableWindow(bMolName);
	this->mButton_MolCreate.EnableWindow(bMolName);
	this->mEdit_Base.EnableWindow(bMolName);
	this->mEdit_TAngle.EnableWindow(bMolName);
}

void CMolOperateDlg::show_read_mol() {
	this->mEdit_MolFilename.EnableWindow(bReadMolFromFile);
	this->m_MolType1.EnableWindow(bReadMolFromFile);
	this->mButton_FileBrowse.EnableWindow(bReadMolFromFile);
	this->mButton_ReadMol.EnableWindow(bReadMolFromFile);
}

void CMolOperateDlg::show_multi_mol() {
	this->mEdit_MultiMolFname.EnableWindow(bMultiMol);
	this->m_MolType2.EnableWindow(bMultiMol);
	this->mButtom_BrowseMultiMolFname.EnableWindow(bMultiMol);
	this->mButton_ReadMultiMol.EnableWindow(bMultiMol);
}

void CMolOperateDlg::show_run_command() {
	this->mCommand_Box.EnableWindow(bRunCommand);
	this->mButton_RunCommand.EnableWindow(bRunCommand);
}

void CMolOperateDlg::show_dialog() {
	CDataExchange dx(this, FALSE); // show data
	CDataExchange *pDX = &dx;

	DDX_Check(pDX, IDC_CHECK_MOLNAME, bMolName);
	DDX_Check(pDX, IDC_CHECK_READFROMFILE, bReadMolFromFile);
	DDX_Check(pDX, IDC_CHECK_READMULTI_MOL, bMultiMol);
	DDX_Check(pDX, IDC_CHECK_RUNCOMMAND, bRunCommand);

	show_mol_name(); show_read_mol(); show_multi_mol(); show_run_command();
}

void CMolOperateDlg::check_mol_name() {
	CDataExchange dx(this, TRUE); // save data
	CDataExchange *pDX = &dx;
	DDX_Check(pDX, IDC_CHECK_MOLNAME, bMolName);
}

void CMolOperateDlg::check_read_from_file() {
	CDataExchange dx(this, TRUE); // save data
	CDataExchange *pDX = &dx;
	DDX_Check(pDX, IDC_CHECK_READFROMFILE, bReadMolFromFile);
}

void CMolOperateDlg::check_read_multi_mol() {
	CDataExchange dx(this, TRUE); // save data
	CDataExchange *pDX = &dx;
	DDX_Check(pDX, IDC_CHECK_READMULTI_MOL, bMultiMol);
}

void CMolOperateDlg::check_run_command() {
	CDataExchange dx(this, TRUE); // save data
	CDataExchange *pDX = &dx;
	DDX_Check(pDX, IDC_CHECK_RUNCOMMAND, bRunCommand);
}

// CMolOperateDlg message handlers

void CMolOperateDlg::OnBnClickedCheckMolname()
{
	// TODO: Add your control notification handler code here
	check_mol_name();
	if (this->bMolName) {this->bReadMolFromFile = FALSE; this->bMultiMol = FALSE; this->bRunCommand = FALSE;}
	show_dialog();
}

void CMolOperateDlg::OnBnClickedCheckReadfromfile()
{
	// TODO: Add your control notification handler code here
	this->check_read_from_file();
	if (this->bReadMolFromFile) {this->bMolName = FALSE; this->bMultiMol = FALSE; this->bRunCommand = FALSE;}
	show_dialog();
}

void CMolOperateDlg::OnBnClickedCheckReadmultiMol()
{
	// TODO: Add your control notification handler code here
	this->check_read_multi_mol();
	if (this->bMultiMol) {this->bMolName = FALSE; this->bReadMolFromFile = FALSE; this->bRunCommand = FALSE;}
	show_dialog();
}

void CMolOperateDlg::OnBnClickedCheckRuncommand()
{
	// TODO: Add your control notification handler code here
	this->check_run_command();
	if (this->bRunCommand) {this->bMolName = FALSE; this->bMultiMol = FALSE; this->bReadMolFromFile = FALSE;}
	show_dialog();
}

void CMolOperateDlg::OnBnClickedButton2()
{
	// TODO: Add your control notification handler code here
	// File browse
	CFileDialog dlgFile(TRUE, LPCTSTR(CString("*.*"))); // to open file
	char fname[256] = "\0";
	if (dlgFile.DoModal() == IDCANCEL) return;
	this->mMolFilename = dlgFile.GetPathName();
	
	CDataExchange dx(this, FALSE); // show 
	CDataExchange *pDX = &dx;
	DDX_Text(pDX, IDC_EDIT_MOL_FNAME, mMolFilename);
}

void CMolOperateDlg::OnBnClickedButtonCreateMolname()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE); // save data
	char polymer_name[256] = "\0";
	char ta_buffer[256] = "\0";
	if (this->mMolName.GetLength() == 0) {
		show_msg("polymer name is not given"); return;
	}
	if (this->mTAngle.GetLength() == 0) {
		show_msg("torsion angle of each monomer is not given"); return;
	}

	string2char(this->mMolName, polymer_name);
	string2char(this->mTAngle, ta_buffer);
	int nMonomer = nMonomerFromPolymerName(polymer_name);
	if (nMonomer < 2) {
		show_msg("polymer has 2 more monomer at least"); return;
	}
	float *ta = new float[nMonomer - 1];
	char *str = ta_buffer;
	if (sscanf(ta_buffer, "%f", ta) != 1) {
		delete[] ta; ta = NULL;
		sprintf(errmsg, "%s ; does not have torsion angle", ta_buffer); return;
	}
	str = ta_buffer;
	for (int i = 1; i < nMonomer - 1; i++) {
		NextSection(&str, ' ');
		if (str == NULL || sscanf(str, "%f", ta + i) != 1) ta[i] = ta[i-1];
	}

	MMOLECULE *m = new MMOLECULE;
	int nbase = this->nBase;
	if (nbase < 0) nbase = 0;
	//ConstructSingleChainPolymerMonomerChain(MMOLECULE &mm, char *polymer, float *ta, float ta_tail, int &nbase)
	if (!ConstructSingleChainPolymerMonomerChain(*m, polymer_name, ta, ta[nMonomer-2], nbase)) {
		delete[] ta; delete m;
		sprintf(errmsg, "failure to construct polymer with polymer name %s and torsion angle %s", polymer_name, ta_buffer);
		show_msg(errmsg);
		return;
	}
	if (nbase > m->nCluster) nbase = m->nCluster - 1;
	m->base = m->cluster + nbase;
	m->base->disable_parent();
	setup_parent_child_relationship(m->base, true); // report error when free hinge exist
	m->checkMolecule(true);

	if (!check_parent_child_setup(m->cluster, m->nCluster)) {
		delete[] ta; delete m;
		return;
	}
	if (!m->linkPars(atompar_db, errmsg)) {
		delete[] ta; delete m;
		show_msg(errmsg); return;
	}
	m->calTorsionAngle();

	m->setup_cluster_matom(); // require the atomic mass, and the parent-child relationship

	// set mass center to {0, 0, 0}
	VECTOR<3> cm;
	m->SetMassCenter(cm.v, true, true);

	LIST<MMOLECULE> *ml = NULL;
	if (mlist == NULL) {ml = new LIST<MMOLECULE>; ml->p = m; mlist = ml;}
	else ml = mlist->add(m, true);

	//::save_xyz_struct(*m, "mol.dat");

	// tell VMM to show the molecule
	int indx = mlist->indx(ml);
	if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);

	nCurMol = mlist->indx(ml);
	sprintf(errmsg, "set molecule # %d as current operating molecule", nCurMol);
	show_msg(errmsg); strcpy(errmsg, "\0");
}

void CMolOperateDlg::OnBnClickedButtonReadmolFromfile()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE); // save data
	if (this->mMolFilename.GetLength() == 0) {show_msg("no file name is given"); return;}
	char fname[256] = "\0";
	string2char(this->mMolFilename, fname);

	switch (this->m_MolType1.GetCurSel()) {
	case 0: // macromol 
		{
		MMOLECULE *m = new MMOLECULE;
		if (!ConstructMoleculeStruct(m, fname)) {
			delete m; m = NULL;
			sprintf(errmsg, "failure to construct molecule with the parameter in file %s", fname);
			show_msg(errmsg); return;
		}

		LIST<MMOLECULE> *ml = NULL;
		if (mlist == NULL) {mlist = new LIST<MMOLECULE>; mlist->p = m; ml = mlist;}
		else ml = mlist->add(m, true);

		int indx = mlist->indx(ml);
		if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);
		nCurMol = indx;
		}
		break;
	case 1: // coarse-grained macromol
		{
		CG_MMOLECULE *cgm = new CG_MMOLECULE;
		if (!ConstructCGMMStruct(cgm, fname)) {
			delete cgm; cgm = NULL;
			sprintf(errmsg, "failure to construct coarse-grained molecule with the parameter in file %s", fname);
			show_msg(errmsg); return;
		}

		LIST<CG_MMOLECULE> *cgml = NULL;
		if (cgmlist == NULL) {cgmlist = new LIST<CG_MMOLECULE>; cgmlist->p = cgm; cgml = cgmlist;}
		else cgml = cgmlist->add(cgm, true);

		int indx = cgmlist->indx(cgml);
		if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);
		nCurMol = indx;
		}
		break;
	default:
		sprintf(errmsg, "unknown macro-mol type");
		show_msg(errmsg); return;
		break;
	}

	sprintf(errmsg, "set molecule # %d as current operating molecule", nCurMol);
	show_msg(errmsg); strcpy(errmsg, "\0");

	return;
}

void CMolOperateDlg::OnBnClickedButton3()
{
	// TODO: Add your control notification handler code here
	// run command
	UpdateData(TRUE); // save data
	if (this->mCommand.GetLength() == 0) {
		show_msg("No command is given"); return;
	}
	string2char(this->mCommand, cmd_msg);
	char c1[50] = "\0", c2[50] = "\0";
	if (sscanf(cmd_msg, "%s %s", c1, c2) == 1) strcpy(c2, "\0");
	::pCurrentCmd = GetCommand(c1, c2);
	if (::pCurrentCmd == NULL) {sprintf(errmsg, "UNKNOWN COMMAND : %s", cmd_msg); show_msg(errmsg); return;}
	pCurrentCmd->func();
	this->mCommand_Box.InsertString(0, this->mCommand);
	int nItems = this->mCommand_Box.GetCount();
	if (nItems > 100) this->mCommand_Box.DeleteString(nItems - 1);
	return;
}

int CMolOperateDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  Add your specialized creation code here
	show = false;
	// all molecules are removed
	release_list<MMOLECULE>(&mlist, true);
	::atompar_db.release(); LJPARS.release();

	return 0;
}

void CMolOperateDlg::OnClose()
{
	// TODO: Add your message handler code here and/or call default
	show = false;
	release_list<MMOLECULE>(&mlist, true);
	::atompar_db.release(); LJPARS.release();
	nCurMol = -1;
	CDialog::OnClose();
}

void CMolOperateDlg::OnBnClickedButtonClearMonitor()
{
	// TODO: Add your control notification handler code here
	release_list<MMOLECULE>(&mlist, true);
	nCurMol = -1;
	if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_CLEAR_MONITOR, 0, 0);
	::atompar_db.release(); LJPARS.release();
}

void CMolOperateDlg::OnBnClickedButtonBrowseMultiMolFname()
{
	// TODO: Add your control notification handler code here
	CFileDialog dlgFile(TRUE, LPCTSTR(CString("*.*"))); // to open file
	char fname[256] = "\0";
	if (dlgFile.DoModal() == IDCANCEL) return;
	this->mMultiMolFname = dlgFile.GetPathName();
	
	CDataExchange dx(this, FALSE); // show 
	CDataExchange *pDX = &dx;
	DDX_Text(pDX, IDC_EDIT_MULTI_MOL_FNAME, mMultiMolFname);
}

void CMolOperateDlg::OnBnClickedButtonReadMultiMolFromfile()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE); // save data
	if (this->mMultiMolFname.GetLength() == 0) {show_msg("no file name is given"); return;}
	char fname[256] = "\0";
	string2char(this->mMultiMolFname, fname);

	switch (this->m_MolType2.GetCurSel()) {
	case 0: // macro-molecules
		{
		release_list<MMOLECULE>(&mlist, true);
		nCurMol = -1;
		if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_CLEAR_MONITOR, 0, 0);

		if (!ReadMultiMol2List_FromFile<MMOLECULE>(&mlist, (void*)(&ConstructMoleculeStruct), fname, "[MACRO-MOLECULE]")) {
			release_list<MMOLECULE>(&mlist, true);
			show_msg(errmsg); return;
		}

		LIST<MMOLECULE> *ml = NULL;
		int indx = 0;
		ml = mlist;
		while (ml != NULL) {
			if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);
			ml = ml->next; indx++;
		}
		}
		break;
	case 1: // coarse-grained macro-molecules
		{
		release_list<CG_MMOLECULE>(&cgmlist, true);
		nCurMol = -1;
		if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_CLEAR_MONITOR, 0, 0);

		/*
		if (!ReadMultiMol2List(fname)) {
			release_list<CG_MMOLECULE>(&cgmlist, true);
			show_msg(errmsg); return;
		}
		*/

		LIST<CG_MMOLECULE> *cgml = NULL;
		int indx = 0;
		cgml = cgmlist;
		while (cgml != NULL) {
			if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);
			cgml = cgml->next; indx++;
		}
		}
		break;
	default:
		sprintf(errmsg, "unknown molecule's type"); show_msg(errmsg); return;
	}

	nCurMol = 0;
	sprintf(errmsg, "set molecule # %d as current operating molecule", nCurMol);
	show_msg(errmsg); strcpy(errmsg, "\0");

	return;
}
