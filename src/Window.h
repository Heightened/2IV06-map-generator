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
    wxDECLARE_EVENT_TABLE();

	Canvas* mapPreview;
	wxGLContext* glContext;

public:
	GeneratorFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

	void InitializeGL();
};

enum {
    ID_New = 1,
	ID_Open = 2,
	ID_Save = 3,
	ID_Export = 4
};
