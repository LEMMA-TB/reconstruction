static const float convToNs = 0.78125;
static const float velWireProp= 24.4;  // (cm/ns)

static const int MinNHit_Phi = 6;
//static const int MaxNHit_Phi = 20;
static const int MaxNHit_Phi = 16; //original

static const int MinNHit_1SLPhi = 4;
static const int MaxNHit_1SLPhi = 8; //original

static const int MinNHit_Theta = 3;
// static const int MinNHit_Theta = 4; //original
static const int MaxNHit_Theta = 8; //original

static const int MaxNHit = 100;
// static const int MaxNHit = 30;

//static const float T0_shift = 8.;
//static const float T_shift = 12.5;
static const float T_shift = 8.3;
static const float T0_shift = 0.;
//static const float T0_max = 40.;
static const float T0_max = 30.;

// sigmaTimes per fit_t0 (o solo phi o solo theta)
// per il fit globale sigmaTimes e' ancora 3.
static const float sigmaTimes_Phi = 5.;
static const float sigmaTimes_1SL = 5.;
