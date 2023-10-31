class COMMAND {
public:
	char cmd1[50];
	char cmd2[50];
	char format[256];
	char help_msg[1024];
	bool (*func) ();
	COMMAND *next;

	COMMAND() {
		strcpy(cmd1, "\0");
		strcpy(cmd2, "\0");
		strcpy(format, "\0");
		strcpy(help_msg, "\0");
		func = NULL;
		next = NULL;
	};
	void Help_String(char* format_string, char* help_msg_string) {
		strcpy(format, format_string);
		strcpy(help_msg, help_msg_string);
	};
	void delCommand(COMMAND *pCmd) {
		if (pCmd == NULL) return;
		if (pCmd->next != NULL) {
			delCommand(pCmd->next);
			pCmd->next = NULL;
		}
		delete pCmd;
	};
	~COMMAND() {
		delCommand(next);
		next = NULL;
	};
};

void UpperString(char* str);

void DefCommands();
COMMAND* GetCommand(char* c1, char* c2);
COMMAND* AddCommand(char* c1, char* c2, bool (*func)());
void DelCommands();
bool ShowCommand();

bool tweak_mol();
bool tweak();
bool shift();
bool copy_cur_mol();
bool show_cur_mol();
bool delete_cur_mol();
bool set_cur_mol();
bool save_cur_mol();
bool highlight();
bool read_charge_mol2();
bool show_charge_cluster();
bool average_charge_cluster();
bool total_dipole();
bool monomer_dipole();
bool monomer_dipole_distribution();
bool cal_cluster_radius();
