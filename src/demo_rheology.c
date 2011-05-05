

typedef enum { STOKES_IVTYPE_CONSTANT, STOKES_IVTYPE_FK, STOKES_IVTYPE_ARRHENIUS,STOKES_IVTYPE_COMPOSITE } StokesIsotropicFunction;
typedef enum { STOKES_RHSTYPE_CONSTANT_G, STOKES_RHSTYPE_VECTOR_G } StokesRHSFunction;

typedef struct Stokes_A_Isotropic {
	StokesIsotropicFunction rheology_type;
	int index;
	/* params for constant rheology */
	double CNST_eta_0;
	/* params for frank-kam rheology */
	double FK_theta, FK_eta_0;
	/* params for arrhenius rheology */
	double AR_E, AR_eta_0;
	/* params for composite rheology */
	int CMP_n,*CMP_phase_list;
	int CMP_comp_type;

};

struct Stokes_RHS {
	StokesRHSFunction rheology_type;
	int index;
	/* params for constant g */
	double CNST_gravity;
	/* params for variable g */
	double V_gravity_mag;
};


void StokesIsotropicViscosityEvaluate_standard(Stokes_A_Isotropic rlist[],DataBucket db)
{
	ParticleTypeStd *basep;
	ParticleTypeViscous *vp;
	int np;
	
	DataBucketGetField(db,"baseparticle",&basep);
	DataBucketGetField(db,"viscousparticle",&vp);
	for (p=0; p<np; p++) {
		phase            = basep[p];
		
		switch(rlist[phase]->rheology_type){
				case STOKES_IVTYPE_CONSTANT:
					vp[p]->viscosity = rlist[phase]->CNST_eta_0;
				break;
			case STOKES_IVTYPE_FK:
				vp[p]->viscosity = rlist[phase]->FK_theta;
				break;
				
				break;
		}
		
	}
}


void StokesIsotropicViscosityEvaluate_director(Stokes_A_Isotropic rlist[],DataBucket db)
{
	ParticleTypeStd      *basep;
	ParticleTypeViscous  *vp;
	ParticleTypeDirector *dvp;
	int np;
	
	DataBucketGetField(db,"baseparticle",&basep);
	DataBucketGetField(db,"director_variables_particle",&dvp);
	DataBucketGetField(db,"viscousparticle",&vp);
	for (p=0; p<np; p++) {
		phase            = basep[p];
		if (rlist[phase]->rheology_type != STOKES_IVTYPE_DIRECTOR) { exit(0); }

		/* compute any auxillary shit */

		/* update director variables */

		/* compute a viscosity */

	}
}

void StokesIsotropicViscosityEvaluate_composite(Stokes_A_Isotropic rlist[],DataBucket db)
{
	ParticleTypeStd      *basep;
	ParticleTypeViscous  *vp;
	int np;
	
	DataBucketGetField(db,"baseparticle",&basep);
	DataBucketGetField(db,"viscousparticle",&vp);
	for (p=0; p<np; p++) {
		phase            = basep[p];
		if (rlist[phase]->rheology_type != STOKES_IVTYPE_COMPOSITE) { exit(0); }
		
		
		
	}
}


void StokesIsotropicViscosityEvaluate(Stokes_A_Isotropic rlist[],DataBucket db,SwarmUpdateRequires sut[])
{
	i = 0;
	while(sut[i]!=NULL) {
		t = sut[i];
		if (t==STANDARD_VARIABLES) {
			
		}
		else if (t==DIRECTOR_VARIABLES) {
			
		}
		else if (t==COMPOSITE_VARIABLES) {
			
			
		}
		else {
			
		}
		i++;
	}
}


