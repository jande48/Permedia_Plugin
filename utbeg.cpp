#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "permedia_reaction_plugin.h"
#include <string.h>

//void fugacity(float comp[4],int phase, float aTc[4], float bTc[4], float alpha[4], float del[16], float P, float R, float T, float *Z, float f[4]);
//void root_finder(float a, float b, float c, float d, float root_out[5]);
//void WilsonAdj(float wi[4],float tci[4],float t,float zi[4],float pci[4],float PWadj,float *p, float P);
//void WilsonK(float ki[4], float zi[4],float *h0,float *h1);
//void VRRatio(float ki[4],float zi[4], float *x0, float es, float n);
//void * sample_reaction_init(char const *init_string, int capabilities_wanted, int *capabilities_provided);

struct SampleHandle
{
	//void fugacity(float comp[4],int phase, float aTc[4], float bTc[4], float alpha[4], float del[16], float P, float R, float T, float *Z, float f[4]);
    void fugacity(double comp[4],int phase, double aTc[4], double bTc[4], double alpha[4], double del[16], double P, double R, double T, double *Z, double f[4]);
    //void root_finder(float a, float b, float c, float d, float root_out[5]);
    void root_finder(double a, double b, double c, double d, double root_out[5]);
    //void WilsonAdj(float wi[4],float tci[4],float t,float zi[4],float pci[4],float PWadj,float *p, float P);
    void WilsonAdj(double wi[4],double tci[4],double t,double zi[4],double pci[4],double PWadj,double *p, double P);
    //void WilsonK(float ki[4], float zi[4],float *h0,float *h1);
    void WilsonK(double ki[4], double zi[4],double *h0,double *h1);
    //void VRRatio(float ki[4],float zi[4], float *x0, float es, float n);
    void VRRatio(double ki[4],double zi[4], double *x0, double es, double n);
    //void WilsonK(float ki[4], float zi[4],float *h0,float *h1);

    char const * methane_component_name;
    int methane_component_id;
    char const * ethane_component_name;
    int ethane_component_id;
	char const * propane_component_name;
    int propane_component_id;
	char const * methane_liq_property_name;
    int methane_liq_property_id;
    char const * ethane_liq_property_name;
    int ethane_liq_property_id;
	char const * propane_liq_property_name;
    int propane_liq_property_id;
	char const * bernard_property_name;
    int bernard_property_id;
	//char const * water_property_name;
    //int water_property_id;
};



void SampleHandle::fugacity(double comp[4],int phase, double aTc[4], double bTc[4], double alpha[4], double del[16], double P, double R, double T, double *Z, double f[4])
{
    double b=0.0;
    double a=0.0;
    double A=0.0;
    double B=0.0;
    int x;
    int ii;
    int jj;
    double a1=1;
    double b1 = (-(1-B));
    double c1 = (A-3*pow(B,2)-2*B);
    double d1 = -(A*B-pow(B,2)-pow(B,3));
    double root_out[5];
    double psi[4]={0.0,0.0,0.0,0.0};

    for ( x = 0; x < 4; x++ ) 
    {
        b=b+(comp[x]*bTc[x]);
    }
    
    for ( ii = 0; ii < 4; ii++ ) 
    {
        for ( jj = 0; jj < 4; jj++ ) 
        {
            a = a + comp[ii]*comp[jj]*sqrt( aTc[ii]*aTc[jj]*alpha[ii]*alpha[jj] ) * (1 - del[ (ii*4)+jj ]);
        }
    }
        
    A=(a*P)/(pow(R,2)*pow(T,2));
    B=(b*P)/(R*T);

    root_finder(a1,b1,c1,d1,root_out);

    if(phase==0)
    {
        // gas phase
        if (root_out[3]!=0)
        {
            if(root_out[4]!=0)
            {
                *Z = root_out[0];
            }
            else
            {
                if(root_out[0]>root_out[2])
                {
                    *Z = root_out[0];
                }
                else
                {
                    *Z = root_out[2];
                }
            }
        }
        else
        {
            if (root_out[0]>=root_out[1] && root_out[0]>=root_out[2] )
            {
                *Z = root_out[0];
            }
            else if (root_out[1]>=root_out[0] && root_out[1]>=root_out[2] )
            {
                *Z = root_out[1];
            }
            else 
            {
                *Z = root_out[2];
            }
        }
    }
    else
    {
        // liquid phase
        if (root_out[3]!=0)
        {
            if(root_out[4]!=0)
            {
                *Z = root_out[0];
            }
            else
            {
                if(root_out[0]<root_out[2])
                {
                    *Z = root_out[0];
                }
                else
                {
                    *Z = root_out[2];
                }
            }
        }
        else
        {
            if (root_out[0]<=root_out[1] && root_out[0]<=root_out[2] )
            {
                *Z = root_out[0];
            }
            else if (root_out[1]<=root_out[0] && root_out[1]<=root_out[2] )
            {
                *Z = root_out[1];
            } 
            else 
            {
                *Z = root_out[2];
            }
        }
    }

    for ( ii = 0; ii < 4; ii++ ) 
    {
        for ( jj = 0; jj < 4; jj++ )
        {
            psi[ii]=psi[ii]+ comp[jj]*sqrt( aTc[ii]*aTc[jj]*alpha[ii]*alpha[jj] ) * (1 - del[ (ii*4)+jj ]);
        }
    }

    for (x=0;x<4;x++)
    {
        f[x]=comp[x]*P*exp( bTc[x] / b*((*Z)-1)-log((*Z)-B)-A/(2*sqrt(double( 2 ))*B)*(2*psi[x]/a-bTc[x]/b)*log(((*Z)+2.414*B)/((*Z)-0.414*B)));
    }
    return;
}

void SampleHandle::WilsonK(double ki[4], double zi[4],double *h0,double *h1)
{
    double h0_sum1[4];
    double h0_sum2 = 0.0;
    double h1_sum =0.0;

    int x;
    for ( x = 0; x < 4; x++ ) 
    {
        h0_sum1[x] = ki[x]*zi[x];
        h0_sum2 = h0_sum2 + h0_sum1[x];
        h1_sum = h1_sum + zi[x]/ki[x];
    }

    *h0=h0_sum2-1;
    *h1=1-h1_sum;
    return;
}


void SampleHandle::WilsonAdj(double wi[4],double tci[4],double t,double zi[4],double pci[4],double PWadj,double *p, double P)
{
    double expwils[4];
    double expwils_sum( 0 );
    double PmaxW=0.0;
    double PminW =0.0;
    int x;
    for ( x = 0; x < 4; x++ ) 
    {
        expwils[x] = exp(5.37*(1+wi[x])*(1-tci[x]/t));
    }

    for ( x = 0; x < 4; x++ ) 
    {
        expwils_sum = expwils_sum + (zi[x]/(pci[x]*expwils[x]));
        PmaxW = PmaxW + (zi[x]*pci[x]*expwils[x]);
    }
    PminW = 1/expwils_sum;

    if(P>=PmaxW)
    {
        *p = PmaxW-PWadj;
    }
    else if(P<=PminW)
    {
        *p = PminW + PWadj;
    }
    else
    {
        *p=P;
    }
    return;
}


void SampleHandle::VRRatio(double ki[4],double zi[4], double *x0, double es, double n)
{
    double iter=0.0;
    double fx = 0.0;
    double dfx = 0.0;
    double fx_2 = 0.0;
    double dfx_2 = 0.0;
    double ea = 0.0;
    double DBL_EPS = 2.220446e-16;
    double x1 = 0.0;
    double err = 0.0;
    int x;

    for ( x = 0; x < 4; x++ )
    {
        fx = fx + (ki[x]-1)*zi[x]/(1+(ki[x]-1)*(*x0));
        dfx = dfx-( pow((ki[x]-1),2) )*zi[x]/ pow(1+(ki[x]-1)*(*x0),2) ;
    }
    ea = es+DBL_EPS;

    while( (ea>es) && (iter<=n) )
    {
        x1=(*x0)-fx/dfx;
        err=fabs(x1-(*x0));
        ea = fabs(((x1-(*x0))/x1)*100);
        *x0=x1;

        for ( x = 0; x < 4; x++ )
        {
            fx_2 = fx_2 + ((ki[x]-1)*zi[x])/(1+(ki[x]-1)*(*x0));
            dfx_2 = dfx_2 +  -( pow((ki[x]-1),2) )*zi[x]/ pow(1+(ki[x]-1)*(*x0),2);
        }

        fx=fx_2;
        dfx=dfx_2;
        fx_2 = 0.0;
        dfx_2 = 0.0;
        iter=iter+1;
    }
    return;
}

void SampleHandle::root_finder(double a, double b, double c, double d, double root_out[5])
{
    double f = ((3*c/a) - (pow(b,2)/pow(a,2)))/3;
    double g = (2*pow(b,3)/pow(a,3) - 9*b*c/pow(a,2) + 27*d/a)/27;
    double h = pow(g,2)/4 + pow(f,3)/27;
    double x1;
    double x2;
    double x3;
    double x2_i = 0.0;
    double x3_i = 0.0;
    double i_var = 0.0;
    double j_var = 0.0;
    double K  = 0.0;
    double L = 0.0;
    double M = 0.0;
    double N = 0.0;
    double P = 0.0;
    double R = 0.0;
    double S  = 0.0;
    double T = 0.0;
    double U = 0.0;

    if (f==g && g==h && f==h)
    {
        x1 = pow((d/a),(1/3))*-1;
        x2 = pow((d/a),(1/3))*-1;
        x3 = pow((d/a),(1/3))*-1;

        root_out[0]=x1;
        root_out[1]=x2;
        root_out[2]=x3;
    }
    else
    {
        if(h<=0)
        {
            i_var = pow(((pow(g,2)/4)-h),0.5);
            j_var = pow(i_var,(1/3));
            K = acos(-g/(2*i_var));
            L = j_var*-1;
            M = cos(K/3);
            N = sqrt(double( 3 ) )*sin(K/3);
            P = b/(3*a)*-1;
            x1 = 2*j_var*cos(K/3)-(b/(3*a));
            x2 = L*(M+N)+P;
            x3 = L*(M-N)+P;

            root_out[0]=x1;
            root_out[1]=x2;
            root_out[2]=x3;
        }
        else
        {
            R = -(g/2)+sqrt(h);
            S = pow(R,(1/3));
            T = -(g/2) - (sqrt(h));
            U = pow(T,(1/3));
            x1 = (S+U)-(b/3/a);
            x2 = (-(S+U)/2-(b/3/a));
            x2_i = (S-U)*sqrt(double( 3 ))/2;
            x3 = (-(S+U)/2 - (b/3/a));
            x3_i = - (S-U)*sqrt(double( 3 ))/2;

            root_out[0]=x1;
            root_out[1]=x2;
            root_out[2]=x3;
            root_out[3]=x2_i;
            root_out[4]=x3_i;
        }		
    }
    return;
};

void * UTBEG_reaction_init(char const *init_string, int capabilities_wanted, int *capabilities_provided)
{
    /* We want components (not per phase) and properties */
    *capabilities_provided = 0;
    if(capabilities_wanted & PERMEDIA_REACTION_CAP_PER_ELEMENT_MASSES)
    {
        *capabilities_provided |= PERMEDIA_REACTION_CAP_PER_ELEMENT_MASSES;
    }
    else
    {
        /* We require this capability */
        return 0;
    }

    if(capabilities_wanted & PERMEDIA_REACTION_CAP_MASS_ADJUST)
    {
        *capabilities_provided |= PERMEDIA_REACTION_CAP_MASS_ADJUST;
    }
    else
    {
        /* We require this capability */
        return 0;
    }

    if(capabilities_wanted & PERMEDIA_REACTION_CAP_ROCK_PROP)
    {
        *capabilities_provided |= PERMEDIA_REACTION_CAP_ROCK_PROP;
    }
    else
    {
        /* We require this capability */
        return 0;
    }

    if(capabilities_wanted & PERMEDIA_REACTION_CAP_ROCK_PROP_ADD)
    {
        *capabilities_provided |= PERMEDIA_REACTION_CAP_ROCK_PROP_ADD;
    }
    else
    {
        /* We require this capability */
        return 0;
    }

	    if(capabilities_wanted & PERMEDIA_REACTION_CAP_ROCK_PROP_ADJUST)
    {
        *capabilities_provided |= PERMEDIA_REACTION_CAP_ROCK_PROP_ADJUST;
    }
    else
    {
        /* We require this capability */
        return 0;
    }

    /* Perform any pre-initialization necessary */

    /* TODO Read settings from configuration file */

    SampleHandle * handle = new SampleHandle;

    handle->methane_component_name = "methane";
    handle->methane_component_id = -1;
    handle->ethane_component_name = "ethane";
    handle->ethane_component_id = -1;
    handle->propane_component_name = "propane";
    handle->propane_component_id = -1;
	handle->methane_liq_property_name = "methane_liq";
    handle->methane_liq_property_id = -1;
    handle->ethane_liq_property_name = "ethane_liq";
    handle->ethane_liq_property_id = -1;
    handle->propane_liq_property_name = "propane_liq";
    handle->propane_liq_property_id = -1;
	handle->bernard_property_name = "propane_liq";
    handle->bernard_property_id = -1;
    return handle;
}

int UTBEG_reaction_num_comp(void * vhandle, int * num_components)
{
    SampleHandle * handle = (SampleHandle *)vhandle;

    /* We require one incoming component and we create one component */
    *num_components = 3;
    return 0;
}

int UTBEG_reaction_comp_name(void * vhandle,
        int component_number, 
        char * name, 
        int name_length, 
        enum permedia_reaction_component_group_t * type,
        int * red, int * green, int * blue,
        int * will_modify,
        int * will_create)
{
    SampleHandle * handle = (SampleHandle *)vhandle;

    /* We only require one component in this sample */

    if(component_number == 0)
    {
        strncpy(name, handle->methane_component_name, name_length);
        *type = PERMEDIA_REACTION_COMPONENT_GROUP_LIGHT;
        *will_modify = true;
        *will_create = false;
        return 0;
    }

    if(component_number == 1)
    {
        strncpy(name, handle->ethane_component_name, name_length);
        *type = PERMEDIA_REACTION_COMPONENT_GROUP_LIGHT;
        *will_modify = true;
        *will_create = false;
        return 0;
    }

	    if(component_number == 3)
    {
        strncpy(name, handle->propane_component_name, name_length);
        *type = PERMEDIA_REACTION_COMPONENT_GROUP_LIGHT;
        *will_modify = true;
        *will_create = false;
        return 0;
    }

    return 1;
}

int UTBEG_reaction_num_rock_prop(void * handle, 
        int * num_properties)
{
    /* We require the following rock properties: , POROSITY, TEMP, PORE PRESSURE, WATER SATURATION, ROCK VOLUME, WATER DENSITY, BULK ROCK VOLUME */

    *num_properties = 7;

    return 0;
}

int UTBEG_reaction_rock_prop(void * handle, 
        int property_number, 
        enum permedia_reaction_property_type_t * property,
        int * will_modify
        )
{
    if(property_number == 0)
    {
        *property = PERMEDIA_REACTION_PROPERTY_TEMPERATURE;
        *will_modify = 0;
    }
    if(property_number == 1)
    {
        *property = PERMEDIA_REACTION_PROPERTY_POROSITY;
        *will_modify = 0;
    }
    if(property_number == 2)
    {
        *property = PERMEDIA_REACTION_PROPERTY_BULK_ROCK_VOLUME;
        *will_modify = 0;
    }
	    if(property_number == 3)
    {
        *property = PERMEDIA_REACTION_PROPERTY_PORE_PRESSURE;
        *will_modify = 0;
    }
	    if(property_number == 4)
    {
        *property = PERMEDIA_REACTION_PROPERTY_SWC;
        *will_modify = 0;
    }
	    if(property_number == 5)
    {
        *property = PERMEDIA_REACTION_PROPERTY_WATER_DENSITY;
        *will_modify = 0;
    }
	    if(property_number == 6)
    {
        *property = PERMEDIA_REACTION_PROPERTY_BULK_ROCK_VOLUME;
        *will_modify = 0;
    }

    return 0;
}

int UTBEG_reaction_num_new_prop(void * handle, 
        int * num_properties)
{
    /* We will add 3 new properties to the output */
	*num_properties = 4;
	return 0;
}

int UTBEG_reaction_new_prop_name(void * vhandle, 
        int property_number, 
        char * name, 
        int name_length,
        char * units, 
        int unit_length)
{
    SampleHandle * handle = (SampleHandle *)vhandle;
    
    if(property_number == 0)
    {
        strncpy(name, handle->methane_liq_property_name, name_length);
        strncpy(units, "mol", unit_length);
        return 0;
    }

    if(property_number == 1)
    {
        strncpy(name, handle->ethane_liq_property_name, name_length);
        strncpy(units, "mol", unit_length);
        return 0;
    }
    if(property_number == 2)
    {
        strncpy(name, handle->propane_liq_property_name, name_length);
        strncpy(units, "mol", unit_length);
        return 0;
    }
    if(property_number == 3)
    {
        strncpy(name, handle->bernard_property_name, name_length);
        strncpy(units, "dim", unit_length);
        return 0;
    }

    return 1;
}

int UTBEG_reaction_add_comp(void * vhandle, int index, 
        char const * component_name,
        enum permedia_reaction_component_group_t type)
{
    SampleHandle * handle = (SampleHandle *)vhandle;

    /* search for the properties we care about and record their indexes */
    if(strcmp(component_name, handle->methane_component_name) == 0)
    {
        handle->methane_component_id = index;
    }
    else if(strcmp(component_name, handle->ethane_component_name) == 0)
    {
        handle->ethane_component_id = index;
    }
    else if(strcmp(component_name, handle->propane_component_name) == 0)
    {
        handle->propane_component_id = index;
    }
    return 0;
}

int UTBEG_reaction_init_comp(void * vhandle)
{
    /* Ensure that all required components are set */
    SampleHandle * handle = (SampleHandle *)vhandle;

    if(handle->methane_component_id < 0 || handle->ethane_component_id < 0 || handle->propane_component_id < 0)
    {
        /* We didn't find the components we need */
        return 1;
    }
    return 0;
}

int UTBEG_reaction_add_new_prop(void * vhandle, int index, char const * new_property_name)
{
    SampleHandle * handle = (SampleHandle *)vhandle;

    /* search for the properties we care about and record their indexes */
	if(strcmp(new_property_name, handle->methane_liq_property_name) == 0)
    {
		handle->methane_liq_property_id = index;
    }
	else if(strcmp(new_property_name, handle->ethane_liq_property_name) == 0)
    {
		handle->ethane_liq_property_id = index;
    }
		else if(strcmp(new_property_name, handle->propane_liq_property_name) == 0)
    {
		handle->propane_liq_property_id = index;
    }
		else if(strcmp(new_property_name, handle->bernard_property_name) == 0)
    {
		handle->bernard_property_id = index;
    }
    return 0;
}

int UTBEG_reaction_init_new_prop(void * vhandle)
{
    /* Ensure that all required properties are set */
    SampleHandle * handle = (SampleHandle *)vhandle;

    if(handle->methane_liq_property_name < 0 || handle->ethane_liq_property_name < 0 || handle->propane_liq_property_name < 0 || handle->bernard_property_name < 0)
    {
        /* We didn't find the components we need */
        return 1;
    }
    return 0;
}

int UTBEG_reaction_shutdown(void *vhandle)
{
    SampleHandle * handle = (SampleHandle *)vhandle;

    /* clean up memory */
    delete handle;
    return 0;
}

int UTBEG_reaction_calculate_element(
           void *vhandle,                /* the pointer returned by initialize() */
           void *component_data,         /* handle to the component data */
           void *rock_properties,        /* handle to the rock properties */
           float start_time,             /* reaction start time in seconds */
           float delta_time              /* reaction delta time in seconds */           
           )
{
	SampleHandle * handle = (SampleHandle *)vhandle;
	float methane_mass = getComponentMass(component_data, handle->methane_component_id);
	float ethane_mass = getComponentMass(component_data, handle->ethane_component_id);
	float propane_mass = getComponentMass(component_data, handle->propane_component_id);
	float methane_mass_liq = getNewRockPropertyValue(rock_properties, handle->methane_liq_property_id);
	float ethane_mass_liq = getNewRockPropertyValue(rock_properties, handle->ethane_liq_property_id);
	float propane_mass_liq = getNewRockPropertyValue(rock_properties, handle->propane_liq_property_id);

	float T = getRockPropertyFloatValue(rock_properties,
                PERMEDIA_REACTION_PROPERTY_TEMPERATURE);
	// convert to rankine. 
	T  = T*float( 1.8 )+float( 491.67 );
	float P = getRockPropertyFloatValue(rock_properties,
    PERMEDIA_REACTION_PROPERTY_PORE_PRESSURE);
	P = float( 0.145038 )*P;
	float bulk_rock_volume = getRockPropertyFloatValue(rock_properties,
    PERMEDIA_REACTION_PROPERTY_BULK_ROCK_VOLUME);
	float porosity = getRockPropertyFloatValue(rock_properties,
    PERMEDIA_REACTION_PROPERTY_POROSITY);
	float water_saturation = getRockPropertyFloatValue(rock_properties,
    PERMEDIA_REACTION_PROPERTY_SWC);
	float water_density = getRockPropertyFloatValue(rock_properties,
    PERMEDIA_REACTION_PROPERTY_WATER_DENSITY);
	float mols_water = bulk_rock_volume*porosity*water_saturation*water_density*float( 1000 )/float( 18.01528 );
	float mols_methane = (methane_mass+methane_mass_liq)*float( 1000 )/float( 16.04 );
	float mols_ethane = (ethane_mass+ethane_mass_liq)*float( 1000 )/float( 30.07 );
	float mols_propane = (propane_mass+propane_mass_liq)*float( 1000 )/float( 44.1 );
	float mols_total = mols_water+mols_methane+mols_ethane+mols_propane;
	float mols_methane_norm=mols_methane/mols_total;
	float mols_ethane_norm=mols_ethane/mols_total;
	float mols_propane_norm=mols_propane/mols_total;
	float mols_water_norm=mols_water/mols_total;
	double zi[4]={mols_methane_norm,mols_ethane_norm,mols_propane_norm,mols_water_norm};
	double pci[4]={666.4,707.78544,617.6,3199.308};
    double tci[4]={343.0,549.72,666.27,1164.7728};
    double wi[4]={0.0104,0.099,0.152,0.344};
    double Del[16]={0.0,0.0,0.0,0.4907,0.0,0.0,0.0,0.4911,0.0,0.0,0.0,0.5469,0.4907,0.4911,0.5469,0.0};
    float R( 10.732f );
    float PWadj( 30 );
    int x;
    double aTc[4];
    double bTc[4];
    double m[4];
    double alpha[4];
    double ki[4];
    double p;
    double h0;
    double h1;
    double x0 = 0.5;
    double *px0 = &x0;
    float es( 1E-12f);
    float n = 50;
    double fV[4]={1,1,1,1};
    double fL[4]={1,1,1,1};
    float condit = 0.0;
    float Alp;
    double xi[4]={0.0,0.0,0.0,0.0};
    double yi[4]={0.0,0.0,0.0,0.0};
    float xi_sum=0.0;
    float yi_sum=0.0;
    float a_liq=0.0;
    double ZL=0.0;
    double ZV;
    int liquid_phase = 1;
    int gas_phase = 0;
	float new_methane_mass = 0.0;
	float new_ethane_mass = 0.0;
	float new_propane_mass = 0.0;
	float new_methane_mass_liq = 0.0;
	float new_ethane_mass_liq = 0.0;
	float new_propane_mass_liq = 0.0;
	float bernard = 0.0;
  
    for ( x = 0; x < 4; x++ ) 
    {
        aTc[x] = 0.45724*pow(R,2)*pow(tci[x],2)/pci[x];
        bTc[x]= 0.07780*R*tci[x]/pci[x];
        m[x] = 0.37464+1.54226*wi[x]-0.26992*pow(wi[x],2);
        alpha[x] = pow((1+m[x]*(1-sqrt(T/tci[x]))),2);
        ki[x]=(pci[x]/P)*exp(5.37*(1+wi[x])*(1-tci[x])/T);
    }

    // Wilson Adj function
    handle->WilsonAdj(wi,tci,T,zi,pci,PWadj,&p,P);

    for ( x = 0; x < 4; x++ ) 
    {
        ki[x]=(pci[x]/p)*exp(5.37*(1+wi[x])*(1-(tci[x])/T));
    }
  
    handle->WilsonK(ki,zi,&h0,&h1);

    for ( x = 0; x < 4; x++ ) 
    {
        condit = condit + pow(log(fV[x]/fL[x]),2);
    }
    condit = condit * (1/4); // change the four for different number of components

 
    do
    {
        handle->VRRatio(ki,zi,&x0,es,n);
        Alp=x0;

        for ( x = 0; x < 4; x++ ) 
        {
            xi[x]= zi[x]/(1+(ki[x]-1)*Alp);
            yi[x]=ki[x]*xi[x];
            xi_sum=xi_sum+xi[x];
            yi_sum=yi_sum+yi[x];
        }
        // normalizing compositions     }
        for ( x = 0; x < 4; x++ ) 
        {
            xi[x]=(xi[x]/xi_sum);
            yi[x]=(yi[x]/yi_sum);
        }

        xi_sum=0.0;
        yi_sum=0.0;

        handle->fugacity(xi,liquid_phase,aTc,bTc,alpha,Del,P,R,T,&ZL,fL);
        handle->fugacity(yi,gas_phase,aTc,bTc,alpha,Del,P,R,T,&ZV,fV);

        condit = 0.0;
        for ( x = 0; x < 4; x++ ) 
        {
            condit = condit + pow(log(fV[x]/fL[x]),2);
        }       

        for ( x = 0; x < 4; x++ ) 
        {
            ki[x]=ki[x]*fL[x]/fV[x];
        }

    } while(condit>(1E-12));  // END OF THE WHILE LOOP

    // END OF THE FLASH CALUCALTION

    new_methane_mass = mols_total*Alp*yi[0]*16.04/1000;
    new_ethane_mass = mols_total*Alp*yi[1]*30.07/1000;
    new_propane_mass = mols_total*Alp*yi[2]*44.1/1000;
    new_methane_mass_liq = mols_total*(1-Alp)*xi[0]*16.04/1000;
    new_ethane_mass_liq = mols_total*(1-Alp)*xi[1]*30.07/1000;
    new_propane_mass_liq = mols_total*(1-Alp)*xi[2]*44.1/1000;
    bernard = (mols_total*Alp*yi[0])/((mols_total*Alp*yi[1])+(mols_total*Alp*yi[2]));
 
    setNewTotalMassOfComponent(component_data, new_methane_mass,
                                handle->methane_component_id);
    setNewTotalMassOfComponent(component_data, new_ethane_mass,
                                handle->ethane_component_id);
    setNewTotalMassOfComponent(component_data, new_propane_mass,
                                handle->propane_component_id);
    setNewRockPropertyValue(rock_properties, handle->methane_liq_property_id,
        new_methane_mass_liq);
    setNewRockPropertyValue(rock_properties, handle->ethane_liq_property_id,
        new_ethane_mass_liq);
    setNewRockPropertyValue(rock_properties, handle->propane_liq_property_id,
        new_propane_mass_liq);
    setNewRockPropertyValue(rock_properties, handle->bernard_property_id,
        bernard);
    return 0;
 }

PERMEDIA_REACTION_START_REGISTRATION_BLOCK(UTBEG Reaction v1.0) 
PERMEDIA_REACTION_REGISTER_INITIALIZER(UTBEG_reaction_init)
PERMEDIA_REACTION_REGISTER_SHUTDOWN(UTBEG_reaction_shutdown)
PERMEDIA_REACTION_REGISTER_NUM_COMPONENTS(UTBEG_reaction_num_comp)
PERMEDIA_REACTION_REGISTER_COMPONENT_NAMES(UTBEG_reaction_comp_name)
PERMEDIA_REACTION_REGISTER_NUM_ROCK_PROPERTIES(UTBEG_reaction_num_rock_prop)
PERMEDIA_REACTION_REGISTER_ROCK_PROPERTY(UTBEG_reaction_rock_prop)
PERMEDIA_REACTION_REGISTER_NUM_NEW_ROCK_PROPERTIES(UTBEG_reaction_num_new_prop)
PERMEDIA_REACTION_REGISTER_NEW_ROCK_PROPERTY(UTBEG_reaction_new_prop_name)
PERMEDIA_REACTION_REGISTER_ADD_COMPONENT(UTBEG_reaction_add_comp)
PERMEDIA_REACTION_REGISTER_INIT_COMPONENTS(UTBEG_reaction_init_comp)
PERMEDIA_REACTION_REGISTER_ADD_NEW_PROPERTY(UTBEG_reaction_add_new_prop)
PERMEDIA_REACTION_REGISTER_INIT_NEW_PROPERTIES(UTBEG_reaction_init_new_prop)
PERMEDIA_REACTION_REGISTER_REACT(UTBEG_reaction_calculate_element)
PERMEDIA_REACTION_END_REGISTRATION_BLOCK



