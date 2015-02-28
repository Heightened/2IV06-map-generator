#include <wx/wxprec.h>

#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif

#include "Canvas.h"

class GeneratorApp : public wxApp {
public:
	virtual bool OnInit();
};

class GeneratorFrame : public wxFrame {
	void OnNew(wxCommandEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnExportMap(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
	void OnGenerate(wxCommandEvent& event);
    wxDECLARE_EVENT_TABLE();

	Canvas* mapPreview;
	wxGLContext* glContext;

public:
	GeneratorFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

	void InitializeGL();
};

enum {
    ID_New = wxID_HIGHEST + 1,
	ID_Open = wxID_HIGHEST + 2,
	ID_Save = wxID_HIGHEST + 3,
	ID_Export = wxID_HIGHEST + 4,
	ID_BtnGenerate = wxID_HIGHEST + 5
};
