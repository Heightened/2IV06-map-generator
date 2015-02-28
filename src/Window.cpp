//int main(int argc, const char* argv[]) {
//	
//}

#include "Window.h"

#ifdef __WINDOWS__
#define CANVAS ShaderCanvas
#else
#include "SimpleCanvas.h"
#define CANVAS SimpleCanvas
#endif

GeneratorFrame::GeneratorFrame(const wxString& title, const wxPoint& pos, const wxSize& size) : wxFrame(NULL, wxID_ANY, title, pos, size) {
	
	//Initialize Menu
    wxMenu* menuFile = new wxMenu;
    menuFile->Append(ID_New, "&New\tCtrl-N", "Reset parameters to generate a new type of map");
	menuFile->Append(ID_Open, "&Open...\tCtrl-O", "Load previously saved parameters");
	menuFile->Append(ID_Save, "&Save...\tCtrl-S", "Save parameters to generate similar maps at a later time");

	menuFile->AppendSeparator();
	menuFile->Append(ID_Export, "&Export...\tCtrl-E", "Export the currently generated map");

    menuFile->AppendSeparator();
    menuFile->Append(wxID_EXIT);

    wxMenu* menuHelp = new wxMenu;
    menuHelp->Append(wxID_ABOUT);

    wxMenuBar* menuBar = new wxMenuBar;
    menuBar->Append( menuFile, "&File" );
    menuBar->Append( menuHelp, "&Help" );

    SetMenuBar(menuBar);

	//Initialize Status bar

    CreateStatusBar();
    SetStatusText("Welcome to Build-an-Isle!");

	//Initialize Map preview

	wxLogStream(NULL);

	std::cout << "let's do this" << std::endl;

	mapPreview = new CANVAS(this, wxSize(this->GetClientSize().GetWidth()-120, this->GetClientSize().GetHeight()));

	//Initialize Generation toolbox

	wxPanel* toolboxPanel = new wxPanel(this, -1, wxPoint(this->GetClientSize().GetWidth()-120,0), wxSize(120, this->GetClientSize().GetHeight()), 0, "Generation Toolbox");

	wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
	sizer->Add(new wxButton(toolboxPanel, ID_BtnGenerate, "Generate Map", wxDefaultPosition, wxSize(120,30), wxSHAPED, wxDefaultValidator, "generateButton"), 0, 0, 0);
	sizer->Add(new wxButton(toolboxPanel, -1, "Export Map", wxDefaultPosition, wxSize(120,30), wxSHAPED, wxDefaultValidator, "exportButton"), 0, 0, 0);
	SetSizer(sizer);
}

void GeneratorFrame::InitializeGL() {
	glContext = mapPreview->GetContext();
}

void GeneratorFrame::OnExit(wxCommandEvent& event) {
    Close(true);
}

void GeneratorFrame::OnAbout(wxCommandEvent& event) {
    wxMessageBox("This program is an island terrain generator expanding on a concept published by Amit Patel.", "About Build-an-Isle", wxOK | wxICON_INFORMATION);
}

void GeneratorFrame::OnNew(wxCommandEvent& event) {
    wxLogMessage("TODO");
}

void GeneratorFrame::OnOpen(wxCommandEvent& event) {
    wxLogMessage("TODO");
}

void GeneratorFrame::OnSave(wxCommandEvent& event) {
    wxLogMessage("TODO");
}

void GeneratorFrame::OnExportMap(wxCommandEvent& event) {
    wxLogMessage("TODO");
}

void GeneratorFrame::OnGenerate(wxCommandEvent& event) {
    mapPreview->GenerateGeometry();
	mapPreview->Refresh(false);
}

wxBEGIN_EVENT_TABLE(GeneratorFrame, wxFrame)
    EVT_MENU(ID_New,   GeneratorFrame::OnNew)
	EVT_MENU(ID_Open,   GeneratorFrame::OnOpen)
	EVT_MENU(ID_Save,   GeneratorFrame::OnSave)
	EVT_MENU(ID_Export,   GeneratorFrame::OnExportMap)
    EVT_MENU(wxID_EXIT,  GeneratorFrame::OnExit)
    EVT_MENU(wxID_ABOUT, GeneratorFrame::OnAbout)
	EVT_BUTTON(ID_BtnGenerate, GeneratorFrame::OnGenerate)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(GeneratorApp);

bool GeneratorApp::OnInit() {
    GeneratorFrame* frame = new GeneratorFrame("Build-an-Isle", wxPoint(100, 100), wxSize(800, 600));
	bool init = frame->Show(true);
	frame->InitializeGL();
	return init;
}
