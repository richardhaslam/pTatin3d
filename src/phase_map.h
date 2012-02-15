
#ifndef __PHASE_MAP_H__
#define __PHASE_MAP_H__

typedef enum { PHASE_MAP_POINT_OUTSIDE=-1, PHASE_MAP_POINT_INSIDE=1 } PhaseMapLocationIndicator;

typedef struct _p_PhaseMap *PhaseMap;
struct _p_PhaseMap {
	double x0,y0,x1,y1;
	double dx,dy;
	int mx,my;
	int nphases;
	int *data;
};

void PhaseMapCreate(PhaseMap *map);
void PhaseMapDestroy(PhaseMap *map);
void PhaseMapGetIndex(PhaseMap pm,const int i,const int j, int *index);
void PhaseMapLoadFromFile(const char filename[],PhaseMap *map);
void PhaseMapGetPhaseIndex(PhaseMap phasemap,double xp[],int *phase);
void PhaseMapCheckValidity(PhaseMap phasemap,int phase,int *is_valid);
void PhaseMapGetMaxPhases(PhaseMap phasemap,int *maxphase);
void PhaseMapViewGnuplot(const char filename[],PhaseMap phasemap);

#endif