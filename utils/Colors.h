


/*
int DY_c = kRed+1;

int ttbar_c = kBlue;
int wt_c = kCyan - 6;
int diboson_c = kGreen+3;
int qcd_c = kRed-7;
int tautau_c = kMagenta;
int gamgam_c = kOrange;
*/



Color_t light_blue = 2088;
Color_t navy_c = 2022;

Color_t DY_c = 2011;
Color_t ttbar_c = navy_c;
Color_t wt_c = 2033;
Color_t diboson_c = 2044;
Color_t qcd_c = 2055;
Color_t tautau_c = 2066;
Color_t gamgam_c = 2077;


// https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
//
TColor *light_blue_co = new TColor(light_blue,  123/255., 202./255., 255./255.);
TColor *navy_co = new TColor(navy_c,  0., 41/225., 186/255.);
TColor *DY_co = new TColor(DY_c, 213./255.,94./255.,0.);
TColor *ttbar_co = navy_co;
TColor *wt_co = new TColor(wt_c,  86./255., 180./255., 233./255.);
TColor *diboson_co = new TColor(diboson_c,  0., 158./255., 115./255.);
TColor *qcd_co = new TColor(qcd_c,  243./255., 168./255., 87./255.);
TColor *tautau_co = new TColor(tautau_c,  19./255., 58./255., 54./255.);
TColor *gamgam_co = new TColor(gamgam_c,  240./255., 228./255., 66./255.);


int DY_style = 1001;
int ttbar_style = 1001;
int wt_style = 1001;
int diboson_style = 1001;
int qcd_style = 1001;
int tautau_style = 1001;
int gamgam_style = 1001;
