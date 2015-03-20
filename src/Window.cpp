//int main(int argc, const char* argv[]) {
//	
//}

#include "Window.h"
#include <wx/valnum.h>

#include "MapIO.h"

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

	gen = new Generator(10000, 10000, 20000);
	mapPreview = new CANVAS(this, wxSize(this->GetClientSize().GetWidth()-120, this->GetClientSize().GetHeight()), gen);

	//Initialize Generation toolbox

	wxPanel* toolboxPanel = new wxPanel(this, -1, wxPoint(this->GetClientSize().GetWidth()-120,0), wxSize(120, this->GetClientSize().GetHeight()), 0, "Generation Toolbox");

	wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);

	//Map shaper
	sizer->Add(new wxStaticText(toolboxPanel, -1, "Map shaper:", wxDefaultPosition, wxDefaultSize));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioMapShaperRadial, "Radial", wxDefaultPosition, wxDefaultSize, wxRB_GROUP));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioMapShaperSquare, "Square"));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioMapShaperRound, "Round"));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioMapShaperCrescent, "Crescent"));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioMapShaperBlob, "RedBlobGames"));

	//Point selector
	sizer->Add(new wxStaticText(toolboxPanel, -1, "Point selector:", wxDefaultPosition, wxDefaultSize));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioPointSelectorRandom, "Random", wxDefaultPosition, wxDefaultSize, wxRB_GROUP));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioPointSelectorHex, "Hexagonal"));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioPointSelectorSquare, "Square"));
	sizer->Add(new wxRadioButton(toolboxPanel, ID_RadioPointSelectorPoisson, "Poisson"));

	//Springs
	sizer->Add(new wxStaticText(toolboxPanel, -1, "Number of springs:", wxDefaultPosition, wxDefaultSize));
	wxIntegerValidator<unsigned long> val(NULL);
	val.SetRange(1, 50001);
	sizer->Add(new wxTextCtrl(toolboxPanel, ID_TextSprings, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, val));

	sizer->Add(new wxButton(toolboxPanel, ID_BtnGenerate, "Generate Map", wxDefaultPosition, wxSize(120,30), wxSHAPED, wxDefaultValidator, "generateButton"), 0, 0, 0);
	sizer->Add(new wxButton(toolboxPanel, ID_BtnExport, "Export Map", wxDefaultPosition, wxSize(120,30), wxSHAPED, wxDefaultValidator, "exportButton"), 0, 0, 0);
	SetSizer(sizer);

	Graph *g = new Graph();
	gen->setPolygonGraph(g);
	Graph *hg = new Graph();
	gen->setHeightGraph(hg);
	gen->start();
}

void GeneratorFrame::InitializeGL() {
	glContext = new wxGLContext(mapPreview, NULL);//mapPreview->GetContext();
	mapPreview->Initialize(glContext);
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
	wxFileDialog* ImportDialog = new wxFileDialog(this, _("Open map"), wxEmptyString, wxEmptyString, _("Map files|*.map"), wxFD_OPEN, wxDefaultPosition);

	// Creates a "open file" dialog with 4 file types
	if (ImportDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		wxString CurrentDocPath = ImportDialog->GetPath();
		FILE * file= fopen(CurrentDocPath.c_str().AsChar(), "rb");
		if (file != NULL) {
			Map::Info* info;
			mapPreview->GenerateGeometry(IO::importMap(file, info));
			mapPreview->Refresh(false);
		} else {
			wxLogMessage("Failed to open file");
		}
	}

	// Clean up after ourselves
	ImportDialog->Destroy();
}

void GeneratorFrame::OnSave(wxCommandEvent& event) {
    wxLogMessage("TODO");
}

void GeneratorFrame::OnExportMap(wxCommandEvent& event) {
	wxFileDialog* ExportDialog = new wxFileDialog(this, _("Export map"), wxEmptyString, wxEmptyString, _("Map files|*.map"), wxFD_SAVE, wxDefaultPosition);

	// Creates a "open file" dialog with 4 file types
	if (ExportDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		wxString CurrentDocPath = ExportDialog->GetPath();
		FILE * file= fopen(CurrentDocPath.c_str().AsChar(), "wb");
		if (file != NULL) {
			IO::exportMap(file, gen->getCenters(), gen->getMapInfo());
		} else {
			wxLogMessage("Failed to open file");
		}
	}

	// Clean up after ourselves
	ExportDialog->Destroy();
}

void GeneratorFrame::OnGenerate(wxCommandEvent& event) {
	gen->reset();
	Graph *g = new Graph();
	gen->setPolygonGraph(g);
	Graph *hg = new Graph();
	gen->setHeightGraph(hg);
    gen->start();
	mapPreview->GenerateGeometry();
	mapPreview->Refresh(false);
}

void GeneratorFrame::OnPointRandom(wxCommandEvent& event) {
	gen->setPointSelectorType(POINTSELECTOR_RANDOM);
}

void GeneratorFrame::OnPointHex(wxCommandEvent& event) {
	gen->setPointSelectorType(POINTSELECTOR_HEX);
}

void GeneratorFrame::OnPointSquare(wxCommandEvent& event) {
	gen->setPointSelectorType(POINTSELECTOR_SQUARE);
}

void GeneratorFrame::OnPointPoisson(wxCommandEvent& event) {
	gen->setPointSelectorType(POINTSELECTOR_POISSON);
}

void GeneratorFrame::OnShapeRadial(wxCommandEvent& event) {
	gen->setMapShaperType(MAPSHAPER_RADIAL);
}

void GeneratorFrame::OnShapeSquare(wxCommandEvent& event) {
	gen->setMapShaperType(MAPSHAPER_SQUARE);
}

void GeneratorFrame::OnShapeBlob(wxCommandEvent& event) {
	gen->setMapShaperType(MAPSHAPER_BLOB);
}

void GeneratorFrame::OnShapeRound(wxCommandEvent& event) {
	gen->setMapShaperType(MAPSHAPER_ROUND);
}

void GeneratorFrame::OnShapeCrescent(wxCommandEvent& event) {
	gen->setMapShaperType(MAPSHAPER_CRESCENT);
}

void GeneratorFrame::OnSpringChange(wxCommandEvent& event) {
	gen->setSpringCount(atol(event.GetString().c_str().AsChar()));
}


wxBEGIN_EVENT_TABLE(GeneratorFrame, wxFrame)
    EVT_MENU(ID_New,   GeneratorFrame::OnNew)
	EVT_MENU(ID_Open,   GeneratorFrame::OnOpen)
	EVT_MENU(ID_Save,   GeneratorFrame::OnSave)
	EVT_MENU(ID_Export,   GeneratorFrame::OnExportMap)
    EVT_MENU(wxID_EXIT,  GeneratorFrame::OnExit)
    EVT_MENU(wxID_ABOUT, GeneratorFrame::OnAbout)
	EVT_BUTTON(ID_BtnGenerate, GeneratorFrame::OnGenerate)
	EVT_BUTTON(ID_BtnExport, GeneratorFrame::OnExportMap)
	EVT_RADIOBUTTON(ID_RadioPointSelectorRandom, GeneratorFrame::OnPointRandom)
	EVT_RADIOBUTTON(ID_RadioPointSelectorHex, GeneratorFrame::OnPointHex)
	EVT_RADIOBUTTON(ID_RadioPointSelectorSquare, GeneratorFrame::OnPointSquare)
	EVT_RADIOBUTTON(ID_RadioPointSelectorPoisson, GeneratorFrame::OnPointPoisson)
	EVT_RADIOBUTTON(ID_RadioMapShaperRadial, GeneratorFrame::OnShapeRadial)
	EVT_RADIOBUTTON(ID_RadioMapShaperSquare, GeneratorFrame::OnShapeSquare)
	EVT_RADIOBUTTON(ID_RadioMapShaperBlob, GeneratorFrame::OnShapeBlob)
	EVT_RADIOBUTTON(ID_RadioMapShaperRound, GeneratorFrame::OnShapeRound)
	EVT_RADIOBUTTON(ID_RadioMapShaperCrescent, GeneratorFrame::OnShapeCrescent)
	EVT_TEXT(ID_TextSprings, GeneratorFrame::OnSpringChange)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(GeneratorApp);

bool GeneratorApp::OnInit() {
    GeneratorFrame* frame = new GeneratorFrame("Build-an-Isle", wxPoint(100, 5), wxSize(1200, 1000));
	bool init = frame->Show(true);
	frame->InitializeGL();
	return init;
}
