#include <stddef.h>
#include <Globals.h>

// -------------------------------------------------------------------------
// global variables
// -------------------------------------------------------------------------
int          psExpertMode_=0;
int          psIOExpertMode_=0;
int          psRSExpertMode_=0;
int          psSamExpertMode_=0;
int          psAnaExpertMode_=0;
int          psOptExpertMode_=0;
int          psPlotTool_=0;
int          psRandomSeed_=-1;
int          psAnalysisInteractive_=0;
int          psFAMaxDataPts_=4000;
int          psConstraintSetOp_=0;
char         *psConfigFileName_=NULL;
PsuadeConfig *psConfig_=NULL;
CommManager  *psCommMgr_=NULL;
