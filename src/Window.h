#include <wx/wxprec.h>

#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif

class GeneratorApp: public wxApp {
public:
	virtual bool OnInit();
};

class GeneratorFrame: public wxFrame {
	void OnNew(wxCommandEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnExportMap(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    wxDECLARE_EVENT_TABLE();
public:
	GeneratorFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
};

enum {
    ID_New = 1,
	ID_Open = 2,
	ID_Save = 3,
	ID_Export = 4
};

wxBEGIN_EVENT_TABLE(GeneratorFrame, wxFrame)
    EVT_MENU(ID_New,   GeneratorFrame::OnNew)
	EVT_MENU(ID_Open,   GeneratorFrame::OnOpen)
	EVT_MENU(ID_Save,   GeneratorFrame::OnSave)
	EVT_MENU(ID_Export,   GeneratorFrame::OnExportMap)
    EVT_MENU(wxID_EXIT,  GeneratorFrame::OnExit)
    EVT_MENU(wxID_ABOUT, GeneratorFrame::OnAbout)
wxEND_EVENT_TABLE()
