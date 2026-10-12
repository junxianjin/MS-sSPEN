/*
 * Agilent Technologies All Rights Reserved.
 * This software contains proprietary and confidential
 * information of Agilent Technologies and its contributors.
 * Use, disclosure and reproduction is prohibited without
 * prior consent.
 */

/********************************************************************** Multi-slice Segmented Excitation sequence **********************************************************************/

#include <standard.h>
#ifndef DPS
#include <sys/time.h>
#endif
#include "sgl.c"

/* Additional GRADIENT/RF structures for Arterial Spin Labelling */
SLICE_SELECT_GRADIENT_T sstag_grad,q2_grad,ps_grad,ss3_grad;
GENERIC_GRADIENT_T q2crush_grad,vs_grad,pscrush_grad;
RF_PULSE_T q2_rf,ps_rf,sge_rf;

/********************************************************************** segment selected sinc_grad and RF **********************************************************************/
SLICE_SELECT_GRADIENT_T region_grad;
REFOCUS_GRADIENT_T      regionr_grad;
RF_PULSE_T region_rf;

/********************************************************************** slice selected sinc_grad and RF **********************************************************************/
SLICE_SELECT_GRADIENT_T slice_grad;
REFOCUS_GRADIENT_T      slicer_grad;
RF_PULSE_T slice_rf;

/* Additional GRADIENT structures for Standard Gradient Echo */
READOUT_GRADIENT_T sgero_grad;
PHASE_ENCODE_GRADIENT_T sgepe_grad;
REFOCUS_GRADIENT_T sgeror_grad;
GENERIC_GRADIENT_T sgespoil_grad;

/* GRADIENT kernel nodes, summary and info */
struct SGL_GRAD_NODE_T *epipk=NULL,*epirk=NULL,*sgek=NULL;
struct SGL_KERNEL_SUMMARY_T *epipkSummary=NULL,*epirkSummary=NULL;
SGL_KERNEL_INFO_T epipR,epipP,epipS,epirR,epirP,epirS,sgeR,sgeP,sgeS;

/* External storage for variables that are likely only calculated for ix=1 */
double *epipeamp;
char   **epipename;
double *flipscale;
double epipTime,epirTime,delayToAcq;
double blipm0,perm0;
double sgegTime,sgedToAcq,sgedFromAcq,sgetr_delay;
double taggradamp;

struct timeval tp;
double time1,time2,time3,time4,time5;
int    rtn;

/********************************************************************** initial_shape_chirp_pulse **********************************************************************/
void shape_chirp_rf (RF_PULSE_T *rf, char rfBase[MAX_STR], char rfName[MAX_STR], double pw, double flip,double rof1, double rof2,double bandwith, double ke)
{
   if ((ix > 1) && !sglarray) return;
   initRf(rf);
   strcpy(rf->pulseBase,rfBase);
   strcpy(rf->pulseName,rfName);
   strcpy(rf->rfcoil,rfcoil);
   rf->rfDuration = pw;
   rf->flip = flip;
   rf->flipmult = 1.0;
   rf->rof1 = rof1;
   rf->rof2 = rof2;

   /* Generate the pulse shape if rf->pulseName is appropriate */
   genRf(rf);

   /* If the pulse shape has not been generated proceed as for init_rf */
   switch (rf->type)
   {
     case RF_NULL:
       rf->rfDuration=pw;
       readRfPulse(rf);

       /* Set actual bandwidth according to pulse duration */
       if ((FP_LT(rf->flip,FLIPLIMIT_LOW)) || (FP_GT(rf->flip,FLIPLIMIT_HIGH)))
         /* use excitation bandwidth */
         rf->bandwidth=bandwith*ke;
       else
         /* use inversion bandwidth */
        rf->bandwidth=bandwith*ke;
       break;
   }

   switch (rf->error)
   {
     case ERR_RF_SHAPE_MISSING:
       sgl_abort_message("ERROR: Can not find RF shape '%s.RF'",rfName);
       break;
     case ERR_RF_HEADER_ENTRIES:
       sgl_abort_message("ERROR rf shape '%s': incorrect header information",rfName);
       break;
   }
   if ((rf->header.rfFraction < 0) || (rf->header.rfFraction > 1))
     sgl_abort_message("ERROR rf shape '%s': RF Fraction must be between 0 and 1",rfName);
}

/********************************************************************** calculation_chirp_pulse **********************************************************************/
 double calc_m_chirp(double bandwidth,double rfDuration,double nv)
{
    double nv0;
    nv0=bandwidth*rfDuration*2;
    return nv0/nv;
}

pulsesequence()
{
/*****************************SPEN initial_variable*********************************/
    double p2_bw,pw2,ke,tpwr2,tpwr2f,pw3,p3_bw,tpwr3,tpwr3f;
    double times,nphase0;
    p2_bw=getval("p2_bw");
    p2=getval("pw2");
    ke=getval("ke");
    tpwr2=getval("tpwr2");
    tpwr2f=getval("tpwr2f");
    p3=getval("pw3");
    p3_bw=getval("p3_bw");
    tpwr3=getval("tpwr3");
    tpwr3f=getval("tpwr3f");

    double delay1=getval("delay1");
    double delay2=getval("delay2");
    double delay3=getval("delay3");
    double delay4=getval("delay4");
    double sign1=getval("sign1");
    char    crusher[256]; //flag to denote the directions of crusher gradients for 180 pulse
    double gcrushx,gcrushy;
    getstr("crusher",crusher);
    double sign_ror,sign_per;
    sign_ror=getval("sign_ror");
    sign_per=getval("sign_per");

/***************************** hard180 RF variable declaration *****************************/
    double tpwr902=getval("tpwr902");
    double p902=getval("p902");
	

/***************************** segment selected RF variable declaration *****************************/
    double region_delay1=getval("region_delay1");
    double region_delay2=getval("region_delay2");
    double region_delay3=getval("region_delay3");
    double region=getval("region");
    double regiontof=getval("regiontof");
    double region_tpwr=getval("region_tpwr");
    double region_tpwrf=getval("region_tpwrf");

/***************************** slice selected RF variable declaration *****************************/
    double slice_delay1=getval("slice_delay1");
    double sllice_delay2=getval("slice_delay2");
    double slice_delay3=getval("slicec_delay3");
    double slice_thk=getval("slice_thk");
    double slicetof;
    double slice_tpwr=getval("slice_tpwr");
    double slice_tpwrf=getval("slice_tpwrf");

  /* Internal variable declarations */
    double  freq90[MAXNSLICE],freq180[MAXNSLICE];
    double  dw,te_delay,tr_delay,*tc_delay,tcmax=0.0,tc0=0.0,xtime;
    double  tau1=0.0,tau2=0.0,del1=0.0,del2=0.0,del3=0.0,del4=0.0,thk2fact,minDELTA=0.0,maxDELTA,grora;
    double  episw,*si,delacq,samplem0,rampm0,skipm0,adjm0,rstadj,rsm0adj,tgrid;
    double  lread,*peamp,roamp,roramp,sign=1.0,ropad,readstretch;
    double  segtr,segtr_delay=0.0,*Mz,*Mxy,*flip,csegt1;
    double  volumes,nrefs,nsges,nsgeps;
    char    segtc[MAXSTR],altread[MAXSTR],cseg[MAXSTR],rseg[MAXSTR],name[MAXSTR],platnpmin[MAXSTR];
    int     nread,nphase,nnav,fulletl,maxkzero,ntraces,seg0;
    int     samplenp,rampnp,platnp,minplatnp,dwint,dwmod,oversample,delint;
    int     ncseg=1,autoscale,SCALEFLIP=FALSE,ixone,tstartup;
    int     GE,SE,DW,IR,RS,EPI,TC,ALT,CSEG,RSEG;
    int     shape90,shape180=0;
    int     i,j,k,l,rrr;
/*spiral variable*/
    int me;
	char ropat[256],ropat1[256],ropat2[256],ropat3[256],ropat4[256];
    getstr("ropat",ropat);getstr("ropat1",ropat1);getstr("ropat2",ropat2);getstr("ropat3",ropat3);getstr("ropat4",ropat4);
	char pepat[256],pepat1[256],pepat2[256],pepat3[256],pepat4[256];
    getstr("pepat",pepat);getstr("pepat1",pepat1);getstr("pepat2",pepat2);getstr("pepat3",pepat3);getstr("pepat4",pepat4);
    double gradro,gradpe,gradro1,gradpe1,gradro2,gradpe2,gradro3,gradpe3,gradro4,gradpe4;
    gradro=getval("gradro");gradro1=getval("gradro1");gradro2=getval("gradro2");gradro3=getval("gradro3");gradro4=getval("gradro4");
    gradpe=getval("gradpe");gradpe1=getval("gradpe1");gradpe2=getval("gradpe2");gradpe3=getval("gradpe3");gradpe4=getval("gradpe4");

  /* Diffusion variables */
    double  *roarr,*pearr,*slarr,roval,peval,slval;
    int     nbval,                  /* total number of bvalues*directions */
            nbro,nbpe,nbsl,         /* bvalues*directions along RO, PE, and SL */
            pol=1;                  /* polarity of 2nd diffusion gradient */
    double  *bro,*bpe,*bsl,         /* b-values along RO, PE, SL */
            *brs,*brp,*bsp,         /* the cross-terms */
            *btrace,                /* the trace */
            max_bval=0,             /* and the maximum trace */
            dcrush=0.0,dgss2=0.0,   /* "delta" for crusher and gss2 gradients */
            Dcrush=0.0,Dgss2=0.0;   /* "DELTA" for crusher and gss2 gradients */

  /* Arterial Spin Labelling variables */
    double  tagthk,addtagthk,taggap,tagdist,ti2,aslti1_delay=0.0,aslTime=0.0,addtheta,rfovof,aslfov;
    double  pssTag[MAXSLICE],freqTag[MAXSLICE],freq90rfov[MAXSLICE];
    double  pq2,flipq2,nq2p,q2thk,q2dist,q2Time=0.0,gcrushq2,tcrushq2,pssQ2[MAXSLICE],freqQ2[MAXSLICE];
    double  ti1,ti1min,aslti2_delay=0.0;
    double  gvs,tdeltavs,bvalvs,pps,flipps,nps,addthk,gcrushps,tcrushps,pssPS[MAXSLICE],freqPS[MAXSLICE];
    double  minpss=0.0,maxpss=0.0,centrepss=0.0,psthk,psTime=0.0;
    char    asl[MAXSTR],tagrev[MAXSTR],q2tips[MAXSTR],q2pat[MAXSTR],q2rf[MAXSTR],minti2[MAXSTR],minti1[MAXSTR];
    char    vascsup[MAXSTR],ps[MAXSTR],pspat[MAXSTR],showps[MAXSTR];
    int     ASL,TAG,OneTag,Q2,Q2RF,VS,PS,RFOV;
    int     shapeTag=0,shapeQ2=0,shapePS=0,shape90rfov=0;

  /* Reference standard gradient echo parameters */
    double  sgesw,sgedw,*sgesi,sgete,sgetr,flipsge,sgedacq;
    char    minsgete[MAXSTR],minsgetr[MAXSTR],sgecentric[MAXSTR];
    int     *petab,inv;
    int     SGE,SGEP;

  /* Real-time variables used in this sequence */
    int  vms_slices = v1;   /* Number of slices */
    int  vms_ctr    = v2;   /* Slice loop counter */
    int  vpe_steps  = v3;   /* Actual number of steps through peloop */
    int  vpe_ctr    = v4;   /* PE loop counter */
    int  vpe_mult   = v5;   /* PE multiplier, ranges from -PE/2 to PE/2 */
    int  vss        = v6;   /* Compressed steady-states */
    int  vssc       = v7;   /* Compressed steady-states */
    int  vtrimage   = v8;   /* Counts down from nt, trimage delay when 0 */
    int  vacquire   = v9;   /* Argument for setacqvar, to skip steady state acquires */
    int  vtrigblock = v10;  /* Number of slices per trigger block */
    int  vtraces    = v11;  /* Number of traces per 'fid' */
    int  vtrace_ctr = v12;  /* Trace counter */
    int  vnphase    = v13;  /* Requested number of PE steps */
    int  vprep      = v14;  /* Prepulse variable */
    int  vq2        = v15;  /* Number of Q2TIPS pulses */
    int  vq2_ctr    = v16;  /* Q2TIPS loop counter */
    int  vps        = v17;  /* Number of ASL presat pulses */
    int  vps_ctr    = v18;  /* ASL presat loop counter */
    //int  vcr1        = v19;     // multiplier along RO
    //int  vcr2        = v20;     // multiplier along PE
    //int  vcr3        = v21;     // multiplier along SL
	
    if (ix==1) ixone=TRUE; else ixone=FALSE;

    tstartup=(int)getval("tstartup");

  /* Enable variable rate sampling */
    setacqmode(WACQ|NZ);

  /* Initialize paramaters */
    get_parameters();
    euler_test(); /* increase ramp time by sqrt(3) if oblique */

    nread=(int)getval("nread");      /* Number of readout points */
    nphase=(int)getval("nphase");    /* Number of phase encode points */
    nnav=(int)getval("nnav");        /* Number of navigator echoes */
    episw=getval("episw");           /* EPI spectral width */
    getstr("platnpmin",platnpmin);   /* Flag to set minimum # readout plateau points in ramp sampling */
    minplatnp=(int)getval("minplatnp"); /* Minimum # points the readout plateau must have in ramp sampling */
    rstadj=getval("rstadj");         /* Ramp sampling timing adjustment */
    rsm0adj=getval("rsm0adj");       /* Ramp sampling moment adjustment */
    tgrid=getval("tgrid");           /* Timing grid, default is 50.0 ns */
    ropad=getval("ropad");           /* Pad time for readout gradient when no ramp sampling */
    readstretch=getval("readstretch"); /* Number of pixels to stretch/compress image in readout */
    grora=getval("grora");           /* Tweak of refocus moment to move zipper up/down (shouldn't be required) */
    volumes=getval("volumes");       /* Number of EPI volumes */
    nrefs=getval("nrefs");           /* Number of EPI reference volumes */
    nsges=getval("nsges");           /* Number of standard gradient echoes */
    nsgeps=getval("nsgeps");         /* Number of standard gradient echoes with prep modules */
    getstr("segtc",segtc);           /* Flag for T2* time correction in segmented acquisitions */
    getstr("altread",altread);       /* Alternate read gradient on alternate segments (only nseg odd and >1) */
    getstr("cseg",cseg);             /* Flag for compressed segments w.r.t. imaging slices */
    segtr=getval("segtr");           /* Compressed segment TR */
    autoscale=(int)getval("autoscale"); /* Flag to automatically scale the flips for compressed segments */
    csegt1=getval("csegt1");         /* The T1 to use in automatic scaling of flips for compressed segments */
    getstr("rseg",rseg);             /* Flag to reverse compressed reference segments */

  /* ASL parameters */
    getstr("asl",asl);               /* ASL flag */
    tagthk=getval("tagthk");         /* ASL tag thickness */
    addtagthk=getval("addtagthk");   /* ASL additional tag thickness for FAIR */
    taggap=getval("taggap");         /* ASL tag gap from imaging slice(s) */
    getstr("tagrev",tagrev);         /* Reverse ASL tag position w.r.t. imaging slices */
    OneTag=(int)getval("onetag");    /* One ASL tag for all slices */
    ti2=getval("ti2");               /* ASL tag inversion time */
    getstr("minti2",minti2);         /* Flag for minimum ASL tag inversion time */
    getstr("q2tips",q2tips);         /* Q2TIPS flag */
    getstr("q2pat",q2pat);           /* Q2TIPS pulse patern */
    pq2=getval("pq2");               /* Q2TIPS pulse duration */
    flipq2=getval("flipq2");         /* Q2TIPS pulse flip angle */
    nq2p=getval("nq2p");             /* Number of Q2TIPS saturation pulses */
    q2thk=getval("q2thk");           /* Q2TIPS slice saturation thickness */
    gcrushq2=getval("gcrushq2");     /* Q2TIPS crushing gradient strength */
    tcrushq2=getval("tcrushq2");     /* Q2TIPS crushing gradient duration */
    getstr("q2rf",q2rf);             /* Q2TIPS RF pulse on/off */
    ti1=getval("ti1");               /* Time from ASL tag to Q2TIPS module */
    getstr("minti1",minti1);         /* Flag for minimum time from ASL tag to Q2TIPS module */
    getstr("vascsup",vascsup);       /* Vascular Suppression flag */
    gvs=getval("gvs");               /* Vascular Suppression gradient strength */
    tdeltavs=getval("tdeltavs");     /* Vascular Suppression "delta" */
    getstr("ps",ps);                 /* ASL in plane presaturation flag */
    getstr("pspat",pspat);           /* ASL in plane presat pulse patern */
    pps=getval("pps");               /* ASL in plane presat pulse duration */
    flipps=getval("flipps");         /* ASL in plane presat pulse flip angle */
    nps=getval("nps");               /* Number of ASL in plane presat pulses */
    addthk=getval("addthk");         /* Additional ASL in plane presat slice thickness */
    gcrushps=getval("gcrushps");     /* ASL in plane presat crushing gradient strength */
    tcrushps=getval("tcrushps");     /* ASL in plane presat crushing gradient duration */
    addtheta=getval("addtheta");     /* dtheta=90 rotates FOV in ASL test mode */
    getstr("showps",showps);         /* Show ASL in plane presaturation in ASL test mode */
    rfovof=getval("rfovof");         /* Rotated FOV additional slice offset in ASL test mode */
    aslfov=getval("aslfov");         /* FOV scaling for standard gradient echo in ASL test mode */

  /* Reference standard gradient echo parameters */
    sgesw=getval("sgesw");           /* Standard gradient echo sw */
    sgete=getval("sgete");           /* Standard gradient echo te */
    sgetr=getval("sgetr");           /* Standard gradient echo tr */
    flipsge=getval("flipsge");       /* Standard gradient echo flip */
    getstr("minsgete",minsgete);     /* Flag for minimum standard gradient echo te */
    getstr("minsgetr",minsgetr);     /* Flag for minimum standard gradient echo tr */
    getstr("sgecentric",sgecentric); /* Flag for standard gradient echo centric-out phase encode */

    GE   = (spinecho[0] == 'n');     /* EPI Gradient Echo */
    SE   = (spinecho[0] == 'y');     /* EPI Spin Echo */
    DW   = (diff[0] == 'y');         /* Diffusion */
    IR   = (ir[0] == 'y');           /* Inversion Recovery */
    RS   = (rampsamp[0] == 'y');     /* Ramp sampling */
    TC   = (segtc[0] == 'y');        /* Segement time correction ... */
    if (nseg<2) TC=FALSE;            /* ... only if there is more than 1 segment */
    ALT  = (altread[0] == 'y');      /* Alternate read gradient on alternate segments ... */
    if (nseg<2) ALT=FALSE;           /* ... only if there is more than 1 segment ... */
    if ((int)nseg%2 == 0) ALT=FALSE; /* ... and an odd # segments */
    CSEG = (cseg[0] == 'y');         /* Compressed segments ... */
    if (SE || nseg<2) CSEG=FALSE;    /* ... only if gradient echo and there is more than 1 segment */
    RSEG = (rseg[0] == 'y');         /* Reverse compressed reference segments */
    ASL  = (asl[0] == 'y');          /* Arterial Spin Labelling */
    Q2   = (q2tips[0] == 'y');       /* Q2TIPS */
    VS   = (vascsup[0] == 'y');      /* Vascular Suppression */
    PS   = (ps[0] == 'y');           /* ASL in plane presaturation */
    TAG  = (tag[0] == 'y');          /* ASL tag pulse, RF on='y'/off='n' */
    Q2RF = (q2rf[0] == 'y');         /* ASL Q2TIPS pulse, RF on='y'/off='n' */
    RFOV = (addtheta == 90.0);       /* Rotate FOV in ASL test mode */

    switch ((int)image)
    {
        case -4: /* Reference standard gradient echo scan with preparation modules */
            EPI=FALSE; SGE=FALSE; SGEP=TRUE;
        break;
        case -3: /* Reference standard gradient echo scan */
            EPI=FALSE; SGE=TRUE; SGEP=FALSE;
        break;
        default: /* Normal EPI scan */
            EPI=TRUE; SGE=FALSE; SGEP=FALSE;
        break;
    }

  /* nread must be a multiple of 4 and nphase a multiple of 2 */
    if (nread%4) { nread +=(4-nread%4); putvalue("nread",nread); }
    if (nphase%2) { nphase++; putvalue("nphase",nphase); }

  /* Check nseg, kzero, nnav and set echo train length */
    if (nseg<1) { nseg=1; putvalue("nseg",nseg); }
    if (nseg>nphase) { nseg=nphase; putvalue("nseg",nseg); }
    fulletl=nphase/(int)nseg;
    if (nphase%(int)nseg) fulletl++;
    maxkzero=nphase/2/(int)nseg;
    if (nphase/2%(int)nseg) maxkzero++;
    if (kzero<1) { kzero=1; putvalue("kzero",kzero); }
    if (kzero>maxkzero) { kzero=maxkzero; putvalue("kzero",kzero); }
    if (nnav<1) { nnav=1; putvalue("nnav",nnav); }
    etl=fulletl-maxkzero+kzero;
    putvalue("etl",etl);

  /* Calculate dwell time and spectral width for EPI */
    dwint=(int)floor(2.5e+6/episw);                      /* requested dw in units of 0.4 us */
    if (2.5e+6/episw-(double)dwint >=0.5) dwint++;       /* correct for requested sw roundoff */
    if (dwint<1) dwint=1;                                /* maximum sw is 2.5 MHz (dw = 0.4 us) */
    oversample=1;                                        /* set no oversample as default */
    dwmod=1;                                             /* default dw is multiple of 0.4 us */
    if (pro != 0.0)
    {                                    /* oversample if readout FOV is offset */
        oversample=2;                                      /* set for oversample */
        dwmod *=2;                                         /* for oversample dw must be a multiple of 0.8 us */
    }
    if (dwint>25*oversample) dwmod *=2;                  /* for sw<100kHz we need an additional x2 (why??) */
    while (dwint%dwmod != 0) dwint++;                    /* set dwint to a suitable multiple of 0.4 us */
    episw=2.5e+6/dwint;                                  /* calculate suitable sw */
    dw=dwint*0.4e-6/oversample;                          /* dwell time */
    sw=episw*oversample;                                 /* oversampled sw */
    putvalue("episw",episw);                             /* return value to parameter set */
    putvalue("oversample",oversample);                   /* return value to parameter set */
    putvalue("sw",sw);                                   /* return value to parameter set */

  /* Calculate dwell time and spectral width for standard gradient echo */
    dwint=(int)floor(2.5e+6/sgesw);                      /* requested dw in units of 0.4 us */
    if (2.5e+6/sgesw-(double)dwint >=0.5) dwint++;       /* correct for requested sw roundoff */
    if ((oversample>1) && (dwint%2)) dwint++;            /* if oversampling make dwint a multiple of 0.8 us */
    sgedw=dwint*0.4e-6/oversample;                       /* standard gradient echo dwell time */
    sgesw=2.5e+6/dwint;                                  /* standard gradient echo sweep width */
    if (FP_GT(sgesw,episw)) sgesw=episw;                 /* sweep width can't be greater than for EPI */
    putvalue("sgesw",sgesw);                             /* return value to parameter set */

  /* Calculate number of samples and set defaults for linear sampling */
    samplenp=oversample*nread/2;                         /* # samples in each readout lobe */
    platnp=samplenp;                                     /* # samples on the readout plateau */
    rampnp=0;                                            /* # samples on each readout ramp */
    if ((si=(double *)malloc(samplenp*sizeof(double))) == NULL) nomem(); /* allocate for sampling intervals */

  /* Calculate readout gradient */
    lread=lro*nread/(nread+2*readstretch);               /* allow readout FOV to be scaled */
    init_readout(&ro_grad,"ro",lread,nread,episw);       /* initialize readout gradient structure */
    ro_grad.pad1=ropad; ro_grad.pad2=ropad;              /* pad the readout gradient according to ropad */
    calc_readout(&ro_grad,NOWRITE,"gro","","");          /* calculate gradient */
    for (i=0;i<samplenp;i++) si[i]=dw;                   /* sampling interval for linear sampling */
    delacq=0.5*(ro_grad.duration-samplenp*dw);           /* from start of readout gradient to acquisition */
 /* Allow ramp sampling for EPI */
  /**********************************************************************************************************************ramp sampling***********************************************/
  
  if (RS) {                                          /* ramp sampling */
    samplem0=dw*oversample*ro_grad.roamp;              /* the moment of each sample without oversampling */
    //rampm0=0.5*ro_grad.tramp*ro_grad.roamp;            /* the moment of the ramp */
    
    rampm0=ro_grad.tramp*ro_grad.roamp;
    
    rampnp=(int)floor(rampm0/samplem0);                /* max # non linear samples we can put on ramp */
    if (platnpmin[0] == 'y') minplatnp=2;              /* min # plateau points if min flag is selected */
    /* Shorten plateau to include ramp sampling equally from ramp up and ramp down */
    platnp=nread/2;
    while ((platnp+2*rampnp > nread/2) && (platnp>minplatnp)) {
      platnp -=2;
    }
    putvalue("minplatnp",platnp);                      /* return min (or selected) # plateau points */
    if (platnp+2*rampnp == nread/2) putCmd("platnpmin = 'y'"); /* if using the min # plateau points check the min box */
    if (platnp == nread/2) {
      platnp *=oversample;                             /* if there are no ramp samples set platnp for oversampling */
      rampnp=0;                                        /* there are no ramp samples */
    } else {                                           /* if there are ramp samples calculate the timings */
      init_readout(&ro_grad,"ro",lread,2*platnp,episw); /* reinitialize readout gradient structure */
      calc_readout(&ro_grad,NOWRITE,"gro","","");      /* recalculate gradient */
      platnp *=oversample;                             /* set platnp for oversampling */
      rampnp=samplenp-platnp;                          /* actual # samples on the ramps */
      rampnp/=2;                                       /* actual # samples on a single ramp */
      samplem0 /=oversample;                           /* correct sample moment for oversampling */
      
      //skipm0=0.5*(ro_grad.m0-samplem0*samplenp);       /* readout moment to skip before sampling */
      //if (skipm0<0.0) skipm0=0.0; 
      skipm0=0.0;                                     /* small -ve imprecision for dw = 2us without middle of first sample */
      delacq=sqrt(2*skipm0/ro_grad.slewRate);          /* delay to start of first sample for an ideal readout gradient */
      delacq+=rstadj;                                  /* allow the start time of sampling to be tweaked by rstadj */
      if (delacq<0.0) {                                /* it can't be less than 0.0 */
        rstadj -=delacq; putvalue("rstadj",rstadj);    /* reset rstadj to minimum */
        delacq=0.0;   
      }
      if (delacq>ro_grad.tramp) delacq=ro_grad.tramp;  /* it can't be more than the readout ramp */
      adjm0=(skipm0-delacq*delacq*ro_grad.slewRate/2)/rampnp; /* m0 adjustment to preserve the # ramp points */
      skipm0=delacq*delacq*ro_grad.slewRate/2;         /* adjust readout moment to skip before sampling */
      adjm0 += rsm0adj*samplem0;                       /* allow m0 adjustment of ramp points */
      
      if (adjm0<-samplem0) adjm0=-samplem0;            /* it can't be adjusted down by more than a sample m0 */
     
      for (i=0;i<samplenp/2;i++) {
        //skipm0+=samplem0+adjm0;                        /* include moment to next sample */
        skipm0 = 0.0;
        si[i]=sqrt(2*skipm0/ro_grad.slewRate)-delacq;  /* time between samples */
        delacq+=si[i];                                 /* adjust delay from start to next sample */
        if (delacq>ro_grad.tramp) {                    /* recalculate if we have reached the plateau */
          delacq-=si[i];                               /* reset the delay */
          si[i]=ro_grad.tramp-delacq;                  /* add time to end of ramp */
          si[i]+=(samplem0+adjm0-0.5*ro_grad.tramp*ro_grad.roamp+delacq*delacq*ro_grad.slewRate/2)/ro_grad.roamp;
          break;
        }
      }
      /* Grid the delays */
      if (tgrid == 0.0) tgrid = 50.0;                  /* default timing grid is 50 ns */
      tgrid=1e9/tgrid;                                 /* scalar to put times on integer grid */
      for (i=0;i<samplenp/2;i++) {                     /* granulate timings to tgrid */
        delint=(int)floor(tgrid*si[i]);
        if (tgrid*si[i]-(double)delint >= 0.5) delint++;
        si[i]=delint/tgrid;
        if (si[i]<dw) si[i]=dw;                        /* trap for any -ve delays */
      }
      for (i=0;i<samplenp/2;i++) si[samplenp-1-i]=si[i]; /* the ramp down is the same */
      delacq=0;
      for (i=0;i<samplenp;i++) delacq +=si[i];         /* sum the sampling intervals */
      delacq=0.5*(ro_grad.duration-delacq);            /* figure the delay to start of acquisition */
      if (delacq<0.0) {                                /* -ve delay due to granulation or m0 adjustment of ramp points */
        if (delacq>-1e-6 && si[0]>dw-delacq) {         /* shouldn't be out by more than 1 us for normal granulation */
          si[0]+=delacq; si[samplenp-1]+=delacq;       /* adjust first and last samples */
        } else {
          if (si[0]>dw-2*delacq) {
            si[samplenp-1]+=2*delacq;                  /* adjust just the last sample */
            warn_message("%s warning: Moment tweak of ramp points is severe!",seqfil);
          } else abort_message("%s error: Moment tweak of ramp points is too severe!",seqfil);
        }
        delacq=0.0;                                    /* set delay to zero */
      }
    } 
  }
  
 
  
/* Set sampling intervals */
    for (i=0;i<samplenp;i++) si[i] -=dw;                 /* subtract dw from the sampling intervals */
    putvalue("platnp",platnp/oversample);                /* # samples on the readout plateau */
    putvalue("rampnp",rampnp/oversample);                /* # samples on each readout ramp */
    putvalue("delacq",delacq);                           /* delay from start of readout to first sample */

  /* Calculate np for the acquisition */
    ntraces=nnav+1+(int)etl;
    np=2*samplenp*ntraces;
    putvalue("np",np);
    putvalue("aqtm",aqtm);

  /* EPI echo spacing */
    esp = ro_grad.duration;
    putvalue("esp",esp);

  /* Calculate readout dephase gradient */
    init_readout_refocus(&ror_grad,"ror");               /* readout dephase */
    ro_grad.m0ref*=grora; /* allow tweak of refocus moment to move zipper up/down (shouldn't be required) */
    calc_readout_refocus(&ror_grad,&ro_grad,WRITE,"gror");
    roramp=ror_grad.amp;

  /* RF Power & Bandwidth Calculations */
    shape_rf(&p1_rf,"p1",p1pat,p1,flip1,rof1,rof1);
    calc_rf(&p1_rf,"tpwr1","tpwr1f");
    shape_rf(&sge_rf,"sge",p1pat,p1,flip1,rof1,rof1);
    calc_rf(&sge_rf,"sgetpwr","sgetpwrf");

  /********************************************************************** Segment_selected RF and grad_initial **********************************************************************/

  /* Calculate slice select gradient */
    init_slice(&region_grad,"region",region);                       /* slice select */
    init_slice_refocus(&regionr_grad,"regionr");                 /* slice refocus */
    calc_slice(&region_grad,&p1_rf,WRITE,"gregion");
    calc_slice_refocus(&regionr_grad,&region_grad,WRITE,"gregionr");
/* 180(2) spen */
    shape_chirp_rf(&p3_rf,"p3",p3pat,p3,flip2,rof1,rof1,p3_bw,ke);
    calc_rf(&p3_rf,"","");
    init_slice_butterfly(&ss3_grad,"ss3",lro*10,gcrush,tcrush);
    calc_slice(&ss3_grad,&p3_rf,WRITE,"gss3");

  /********************************************************************** Slice_selected RF and grad_initial **********************************************************************/

  /* Calculate slice select gradient */
    init_slice(&slice_grad,"slice_thk",slice_thk);                       /* slice select */
    init_slice_refocus(&slicer_grad,"slicer");                 /* slice refocus */
    calc_slice(&slice_grad,&p1_rf,WRITE,"gslice");
    calc_slice_refocus(&slicer_grad,&slice_grad,WRITE,"gslicer");

 /********************************************************************** PE_grad_initial **********************************************************************/
    times=calc_m_chirp(p2_bw,p2,nphase);
    nphase0=(int)nphase*times;
  /* Calculate phase encode gradient */
    init_phase(&per_grad,"per",lpe/times,nphase);              /* phase encode offset */
    calc_phase(&per_grad,NOWRITE,"","");                 /* calculate phase encode as usual */

  /* Calculate phase encode blip to be a multiple of 8 usec */
    if (ixone) blipm0=per_grad.m0*2/nphase; /* only do for ix=1 as we may not recalculate per_grad */
    init_dephase(&pe_grad,"pe");
    calc_dephase(&pe_grad,NOWRITE,blipm0*nseg,"","");
    pe_grad.duration=granularity(pe_grad.duration,GRADIENT_RES*2);
    pe_grad.calcFlag=AMPLITUDE_FROM_MOMENT_DURATION;
    calc_dephase(&pe_grad,NOWRITE,blipm0*nseg,"gpe","tpe");
  /* Adjust phase encode offset for kzero */
    if (ixone) per_grad.m0 *= (1-(maxkzero-kzero)*nseg*2/nphase); /* only for ix=1 */
    calc_phase(&per_grad,NOWRITE,"gper","tper");
    per_grad.calcFlag=AMPLITUDE_FROM_MOMENT_DURATION_RAMP;
    if (ixone) perm0 = per_grad.m0; /* only for ix=1 */

  /* Intra-TE gradients */
    if (SE)
    { /* EPI Spin Echo */
        shape_chirp_rf(&p2_rf,"p2",p2pat,p2,flip2,rof1,rof1,p2_bw,ke);
        calc_rf(&p2_rf,"","");
        init_slice_butterfly(&ss2_grad,"ss2",lpe*10,gcrush,tcrush);
        calc_slice(&ss2_grad,&p2_rf,WRITE,"gss2");
        displayRf(&p2_rf);
        displaySlice(&ss2_grad);

        if (strchr(crusher,'x')||strchr(crusher,'y'))
        {
            init_generic(&spoil_grad,"spoil",gcrush,tcrush);
            calc_generic(&spoil_grad,WRITE,"","");
            displayGeneric(&spoil_grad);
            if (strchr(crusher,'x')) gcrushx=1;else gcrushx=0;
            if (strchr(crusher,'y')) gcrushy=1;else gcrushy=0;
        }
    }

  /* Segment time correction */
    tc_delay=0;
    if (TC)
    {
        tc_delay = (double *)malloc(nseg*sizeof(double));  /* Time correction delays */
        for (i=0;i<nseg;i++) tc_delay[i]=granularity(i*ro_grad.duration/nseg,GRADIENT_RES);
        tcmax=tc_delay[(int)nseg-1];                       /* Maximum time correction delay */
        seg0=(nphase/2-1)%(int)nseg;                       /* The segment kzero is in */
        tc0=tc_delay[seg0];                                /* Time correction for the segment kzero is in */
    }

  /* Min TE */
    if (GE)
    { /* EPI gradient echo */}
    else
    {
        tau1 = ss_grad.rfCenterBack+ssr_grad.duration+ss2_grad.rfCenterFront+4e-6;
        if (VS) tau1 += 2*vs_grad.duration;
        tau2 = ss2_grad.rfCenterBack+ror_grad.duration+tc0+(nnav+1+kzero)*esp-esp/2;
        if (DW) {}
        else
        {
            temin = 2*MAX(tau1,tau2);
        }
    }
    if (minte[0] == 'y')
    {
        te = temin;
        putvalue("te",te+1.25e-8);
    }
    if (te < temin)
    {
        abort_message("TE too short, minimum TE= %.2fms\n",temin*1000+0.005);
    }
  /* Calculate delays */
    if (GE) {}
    else
    {
        if (DW) { /* optimally, the 2nd diff grad is right after the 180 */}
        else
        {
            del1 = te/2.0-tau1; del2 = te/2.0-tau2;
        }
        epipTime=0.0;
    }

  /* Put all readout gradients in a kernel */
    if (ixone || sglarray)
    {
        start_kernel(&epirk); /* start epi read kernel */
        add_gradient((void *)&ro_grad,"readNav",READ,START_TIME,"",0.0,PRESERVE);
        for (i=1;i<nnav;i++)
     	add_gradient((void *)&ro_grad,"readNav",READ,BEHIND_LAST,"",0.0,i%2);		
        add_gradient((void *)&ro_grad,"read",READ,BEHIND_LAST,"",0.0,i++%2);
        add_gradient((void *)&per_grad,"phaseDeph",PHASE,SAME_END,"read",0.0,PRESERVE);
        add_gradient((void *)&ro_grad,"read",READ,BEHIND_LAST,"",0.0,i++%2);
		
        for (j=1;j<etl;j++)
        {
            add_gradient((void *)&pe_grad,"phase",PHASE,BEHIND_LAST,"",-pe_grad.duration/2,INVERT);
            add_gradient((void *)&ro_grad,"read",READ,BEHIND_LAST,"",-pe_grad.duration/2,i++%2);
        }
    }

  /* Write the composite gradients */
    peamp  = (double *)malloc(nseg*sizeof(double));        /* PE amplitudes for each shot */
    if (ixone)
    { /* first array element */
        epipeamp  = (double *)malloc(nseg*sizeof(double));   /* external storage of PE amplitudes */
        epipename = (char **)malloc(nseg*sizeof(char *));    /* external storage of PE names */
        for (i=0;i<nseg;i++) epipename[i] = (char *)malloc(20*sizeof(char));
    }
    for (i=0;i<nseg;i++)
    { /* loop over segments */
        if (ixone || sglarray)
        {
            per_grad.m0 = perm0-(i+1)*blipm0;                  /* set PE dephase moment for the shot */
            calc_phase(&per_grad,NOWRITE,"","");               /* calculate PE dephase */
            sprintf(name,"%s%d","epirpe",i+1);                 /* set PE name for the shot */
            if (evaluateSummary(epirk,&epirkSummary))
            {        /* check if the read kernel has been modified */
                set_comp_info(&epirR,"epirro");                  /* set read shape name */
                set_comp_info(&epirP,name);                      /* set PE shape name */
                set_comp_info(&epirS,"epirsl");                  /* set slice shape name */
                epirTime=write_comp_grads(&epirR,&epirP,&epirS); /* write the gradient shapes */
                epipeamp[i]=epirP.amp;                           /* store PE amplitude externally */
                strcpy(epipename[i],epirP.name);                 /* store PE name externally */
            }
        }
        peamp[i]=epipeamp[i]; /* reset local PE amplitude (can get set to 0.0 for reference scans) */
    }
    roamp=epirR.amp; /* reset local read amplitude (can get reversed for reference scan) */

    switch ((int)image)
    {
        case 1: /* Real image scan */
            break;
        case 0: /* Normal reference scan */
            for (i=0;i<nseg;i++) peamp[i]=0.0;
            break;
        case -1: /* Inverted image scan */
            roamp=-epirR.amp; roramp=-ror_grad.amp;
            break;
        case -2: /* Inverted reference scan */
            for (i=0;i<nseg;i++) peamp[i]=0.0;
            roamp=-epirR.amp; roramp=-ror_grad.amp;
            break;
    }

  /* Standard Gradient Echo */
    if ((sgesi=(double *)malloc(samplenp*sizeof(double))) == NULL) nomem(); /* allocate for sampling intervals */
    init_readout(&sgero_grad,"sgero",lro*aslfov,nread,sgesw);   /* initialize readout gradient structure */
    calc_readout(&sgero_grad, NOWRITE,"","","");           /* calculate gradient */
    for (i=0;i<samplenp;i++) sgesi[i]=sgedw-dw;            /* sampling interval for linear sampling */
    sgedacq=0.5*(sgero_grad.duration-samplenp*sgedw);      /* from start of readout gradient to acquisition */
    init_readout_refocus(&sgeror_grad,"sgeror");           /* readout dephase */
    calc_readout_refocus(&sgeror_grad,&sgero_grad,NOWRITE,"");
    init_phase(&sgepe_grad,"sgepe",lpe*aslfov,nphase);     /* phase encode */
    calc_phase(&sgepe_grad,NOWRITE,"","");                 /* calculate phase encode as usual */
    pe_steps=ntraces*nseg;                                 /* actual number of steps through peloop */

  /* Add gradient structures to a kernel */
    if (!EPI || ixone)
    {
        xtime = ss_grad.tramp+ssr_grad.duration;
        xtime = xtime>sgeror_grad.duration ? xtime : sgeror_grad.duration;
        xtime = xtime>sgepe_grad.duration ? xtime : sgepe_grad.duration;
        xtime -= ss_grad.tramp+ssr_grad.duration;
        start_kernel(&sgek); /* start sge kernel */
        add_gradient((void *)&ss_grad,"slice",SLICE,START_TIME,"",0.0,PRESERVE);
        add_gradient((void *)&ssr_grad,"sliceReph",SLICE,BEHIND_LAST,"slice",0.0,INVERT);
        add_gradient((void *)&sgeror_grad,"readDeph",READ,SAME_END,"sliceReph",xtime,INVERT);
        init_dephase(&sgespoil_grad,"sgespoil");
        calc_dephase(&sgespoil_grad,NOWRITE,sgedw*nread*oversample*sgero_grad.amp-sgero_grad.m0def,"","");
        xtime = sgero_grad.tramp;
        xtime = xtime<sgespoil_grad.tramp ? xtime : sgespoil_grad.tramp;
        add_gradient((void *)&sgepe_grad,"phase",PHASE,SAME_END,"readDeph",0.0,PRESERVE);
        add_gradient((void *)&sgero_grad,"read",READ,BEHIND,"readDeph",0.0,PRESERVE);
        add_gradient((void *)&sgepe_grad,"phaseRew",PHASE,BEHIND,"read",0.0,INVERT);
        add_gradient((void*)&sgespoil_grad,"spoil",READ,BEHIND,"read",-xtime,PRESERVE);

        temin = get_timing(FROM_RF_CENTER_OF,"slice",TO_ECHO_OF,"read");
        if (minsgete[0] == 'y')
        {
            sgete = temin;
            putvalue("sgete",sgete+1.25e-8);
        }
        if (sgete < temin)
        {
            abort_message("Standard Gradient Echo TE too short, minimum TE= %.2fms\n",temin*1000+0.005);
        }
    /* Calculate delays */
        te_delay = granularity(sgete-temin,GRADIENT_RES);
        if (te_delay>0.0) change_timing("readDeph",te_delay);
        sgedToAcq = get_timing(FROM_RF_PULSE_OF,"slice",TO_START_OF,"read");
        sgedFromAcq = get_timing(FROM_END_OF,"read",TO_END_OF,"phaseRew");
        xtime = get_timing(FROM_END_OF,"read",TO_END_OF,"spoil");
        sgedFromAcq = sgedFromAcq>xtime ? sgedFromAcq : xtime;
    /* Write the composite gradients */
        set_comp_info(&sgeR,"sgeread");
        set_comp_info(&sgeP,"sgephase");
        set_comp_info(&sgeS,"sgeslice");
        sgegTime=write_comp_grads(&sgeR,&sgeP,&sgeS);
    /* Min TR */
        trmin = sgegTime + 4e-6; /* ensure that sgetr_delay is at least 4us */
        if (!SGEP) trmin *=ns;
        if (minsgetr[0] == 'y')
        {
            sgetr = trmin;
            putvalue("sgetr",sgetr+1.25e-8);
        }
        if (sgetr < trmin)
        {
            abort_message("Standard Gradient Echo TR too short, minimum TR = %.2fms\n",trmin*1000+0.005);
        }
        if (SGE) sgetr_delay = granularity((sgetr-trmin)/ns,GRADIENT_RES);
        else sgetr_delay = granularity(sgetr-trmin,GRADIENT_RES);
    }

  /* Set up frequency offset pulse shape list */
    offsetlist(pss,ss_grad.ssamp,0,freq90,ns,'c');
    shape90 = shapelist(p1_rf.pulseName,ss_grad.rfDuration,freq90,ns,ss_grad.rfFraction,'c');

  /* Calculate spoiler */
    if (tspoil>0.0)
    {
        init_generic(&spoil_grad,"spoil",gspoil,tspoil);
        calc_generic(&spoil_grad,WRITE,"","");
    }

  /* Create optional prepulse events */
    if (sat[0] == 'y')  create_satbands();
    if (fsat[0] == 'y') create_fatsat();
    if (mt[0] == 'y')   create_mtc();
    if (ir[0] == 'y')   create_inversion_recovery();

    VS=0; /* there will be no vascular suppression if we are not doing ASL */
    PS=0; /* there will be no in-plane presaturation if we are not doing ASL */

  /* Calculate B values */
    if (ixone)
    {
  /* Get and check the size of the dro,dpe,dsl arrays */
        P_getVarInfo(CURRENT,"dro",&dvarinfo); nbro = (int)dvarinfo.size; nbval = nbro;
        P_getVarInfo(CURRENT,"dpe",&dvarinfo); nbpe = (int)dvarinfo.size; if (nbpe > nbval) nbval = nbpe;
        P_getVarInfo(CURRENT,"dsl",&dvarinfo); nbsl = (int)dvarinfo.size; if (nbsl > nbval) nbval = nbsl;
        if ((nbro != nbval) && (nbro != 1))
            abort_message("%s: Number of directions/b-values must be the same for all axes (readout)",seqfil);
        if ((nbpe != nbval) && (nbpe != 1))
            abort_message("%s: Number of directions/b-values must be the same for all axes (phase)",seqfil);
        if ((nbsl != nbval) && (nbsl != 1))
            abort_message("%s: Number of directions/b-values must be the same for all axes (slice)",seqfil);
        if ((roarr=(double *)malloc(nbro*sizeof(double))) == NULL) nomem();
        if ((pearr=(double *)malloc(nbpe*sizeof(double))) == NULL) nomem();
        if ((slarr=(double *)malloc(nbsl*sizeof(double))) == NULL) nomem();
        S_getarray("dro",roarr,nbro*sizeof(double)); /* NB S_getarray returns # elements */
        S_getarray("dpe",pearr,nbpe*sizeof(double));
        S_getarray("dsl",slarr,nbsl*sizeof(double));
    /* Allocate for bvalues */
        if ((bro=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((bpe=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((bsl=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((brs=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((brp=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((bsp=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
        if ((btrace=(double *)malloc(nbval*sizeof(double))) == NULL) nomem();
    /* Initialize variables */
        if (SE) { dgss2 = p2/2; Dgss2 = dgss2; dcrush = tcrush; Dcrush = dcrush + p2; }
        roval=roarr[0]; peval=pearr[0]; slval=slarr[0];
    /* Calculate the values */
        for (i = 0; i < nbval; i++)
        {
            if (nbro>1) roval=roarr[i]; if (nbpe>1) peval=pearr[i]; if (nbsl>1) slval=slarr[i];
        /* Readout */
            bro[i]  = bval(gdiff*roval,tdelta,tDELTA);
            bro[i] += (nnav+1+kzero)*bval(ro_grad.roamp,esp/2.0,esp/2.0);
            /* Phase */
            bpe[i]  = bval(gdiff*peval,tdelta,tDELTA);
            /* Slice */
            bsl[i]  = bval(gdiff*slval,tdelta,tDELTA);
            /* Readout/Slice Cross-terms */
            brs[i]  = bval2(gdiff*roval,gdiff*slval,tdelta,tDELTA);
            /* Readout/Phase Cross-terms */
            brp[i]  = bval2(gdiff*roval,gdiff*peval,tdelta,tDELTA);
            /* Slice/Phase Cross-terms */
            bsp[i]  = bval2(gdiff*slval,gdiff*peval,tdelta,tDELTA);
            if (SE)
            {
                /* We must also consider slice gradients & crushers for slice gradients */
                /* Everything else is outside diffusion gradients, and thus constant    */
                /* for all b-values/directions                                          */
                /* Slice */
                bsl[i] += bval(ss2_grad.ssamp,dgss2,Dgss2);
                bsl[i] += bval(gcrush,dcrush,Dcrush);
                bsl[i] += bval_nested(gcrush,dcrush,Dcrush,ss2_grad.ssamp,dgss2,Dgss2);
                bsl[i] += bval_nested(gdiff*slval,tdelta,tDELTA,gcrush,dcrush,Dcrush);
                bsl[i] += bval_nested(gdiff*slval,tdelta,tDELTA,ss2_grad.ssamp,dgss2,Dgss2);
                /* Readout/Slice Cross-terms */
                brs[i] += bval_cross(gdiff*roval,tdelta,tDELTA,gcrush,dcrush,Dcrush);
                brs[i] += bval_cross(gdiff*roval,tdelta,tDELTA,ss2_grad.ssamp,dgss2,Dgss2);
                /* Slice/Phase Cross-terms */
                bsp[i] += bval_cross(gdiff*peval,tdelta,tDELTA,gcrush,dcrush,Dcrush);
                bsp[i] += bval_cross(gdiff*peval,tdelta,tDELTA,ss2_grad.ssamp,dgss2,Dgss2);
            }
            btrace[i] = (bro[i]+bsl[i]+bpe[i]);
            if (max_bval < btrace[i]) max_bval = btrace[i];
        }  /* End for-all-directions */
    /* Write the values and free memory */
        putarray("bvalrr",bro,nbval); putarray("bvalpp",bpe,nbval);
        putarray("bvalss",bsl,nbval); putarray("bvalrp",brp,nbval);
        putarray("bvalrs",brs,nbval); putarray("bvalsp",bsp,nbval);
        putarray("bvalue",btrace,nbval); putvalue("max_bval",max_bval);
        free(bro); free(bpe); free(bsl); free(brs); free(brp); free(bsp);
        free(roarr); free(pearr); free(slarr); free(btrace);
    }

  /* Check nsblock, the number of slices blocked together
     (used for triggering and/or inversion recovery) */
    check_nsblock();

  /* Min TR */
    trmin = epipTime + epirTime + tcmax + tspoil + 12e-6;  /* ensure that tr_delay is at least 4us */
    if (SE || DW || VS) trmin += ss_grad.duration+ssr_grad.duration+del1;
    if (SE) trmin += ss2_grad.duration+del2;
    if (DW) trmin += 2*diff_grad.duration+del3+del4;
    if (VS) trmin += 2*vs_grad.duration;

    if ((sat[0] == 'y') && !ASL) trmin += satTime;
    if ((mt[0] == 'y') && !ASL)  trmin += mtTime;
    if ((fsat[0]=='y') && !ASL)  trmin += fsatTime;
    if (ticks > 0) trmin += 4e-6;

  /* Adjust for all slices */
    trmin *= ns;


    if (mintr[0] == 'y')
    {
        tr = trmin;
        putvalue("tr",tr+1.25e-8);
    }
    if (tr < trmin)
    {
        abort_message("TR too short, minimum TR = %.2fms\n",trmin*1000+0.005);
    }

  /* Calculate tr delay */
    tr_delay = granularity((tr-trmin)/ns,GRADIENT_RES);

  /* For NIFTI slice duration */
    if (trtype) putvalue("niftitslice",trmin/ns);
    else putvalue("niftitslice",trmin/ns+tr_delay);

    sgl_error_check(sglerror);

  /* Return some gradient durations */
    putvalue("rotramp",ro_grad.tramp);
    putvalue("tro",ro_grad.duration);
    putvalue("tror",ror_grad.duration);
    putvalue("petramp",pe_grad.tramp);
    putvalue("tss",ss_grad.duration);
    putvalue("tssr",ssr_grad.duration);
    if (SE) putvalue("tss2",ss2_grad.duration);

  /* Set DDR offset */
    roff=0; /* we oversample if we want pro !=0 */


    if (ixone)
    { /* Attempt to set reasonable experiment time */
        xtime=4e-6;
        if (sat[0] == 'y') xtime += satTime;
        if (mt[0] == 'y')  xtime += mtTime;
        if ((fsat[0]=='y') && !IR && !ASL) xtime += fsatTime;
        if (IR) xtime += irgradTime+4e-6+ti-timin;
    /* Standard gradient echoes with preparation modules */
        xtime = nsgeps*ns*(tr_delay+xtime+(pe_steps+fabs(ssc))*(sgegTime+sgetr_delay+108e-6));
    /* Standard gradient echoes */
        if (ssc<0) xtime += nsges*(ns*tr_delay+ns*(pe_steps-ssc)*(sgegTime+sgetr_delay+108e-6));
        else
        {
            xtime += ns*(tr_delay+(pe_steps+ssc)*(sgegTime+sgetr_delay+108e-6));
            xtime += (nsges-1)*(ns*tr_delay+ns*pe_steps*(sgegTime+sgetr_delay+108e-6));
        }
    /* EPI */
        if (ss<0) xtime += ((nt-ss)*tr*nseg*nrefs+(ntmean-ss)*trmean*nseg*volumes);
        else xtime += (tr*nseg*ss+tr*ntmean*nseg*nrefs+trmean*ntmean*nseg*volumes);
    /* Inter image delay */
        xtime += trimage*(nsgeps+nsges+volumes+nrefs);
        xtime += 1.0; /* add 1 sec for startup */
        g_setExpTime(xtime);
    }

 /********************************************************************** pulse sequence begin **********************************************************************/
  /* PULSE SEQUENCE ****************************************/
    status(A);
    rotate();
    triggerSelect(trigger);
    obsoffset(resto);
    delay(4e-6);
    delay(tr_delay);

  /* Trigger */
    if (ticks > 0) F_initval((double)nsblock,vtrigblock);

  /* Prepulses */
    assign(zero,vprep); /* by default we will always do preparation modules */

  /* Set grad_advance for first array, first acquired transient */
    if (ixone) { ifzero(rtonce); grad_advance(tep); endif(rtonce); }

    for (i=0;i<4;i++)
    { /* loop over standard segments */

        if (EPI && trtype) delay(ns*tr_delay); /* EPI relaxation delay */

    /*************************************************90segment_selected_sinc pulse******************************************************/
        obsoffset(regiontof);
        obspower(p1_rf.powerCoarse);
        obspwrf(p1_rf.powerFine);
        delay(GRADIENT_RES-MIN_EVENT_TIME);
        delay(GRADIENT_RES-MIN_EVENT_TIME);


        obl_shapedgradient(region_grad.name,region_grad.duration,0.0,0.0,region_grad.amp,NOWAIT);
        delay(region_grad.rfDelayFront);

        shaped_pulse(p1pat,region_grad.rfDuration,oph,rof1,rof2);
        delay(region_grad.rfDelayBack);
        obl_shapedgradient(regionr_grad.name,regionr_grad.duration,0,0,-regionr_grad.amp,WAIT);

        obsoffset(tof);
        //delay(region_delay1);

        if (SE)
        {

	    region_delay2=p2;
           // delay(region_delay2);
            delay(delay3);
        if (fsat[0] == 'y') fatsat();

        /*************************************************************180 chirp pulse*****************************************************/
            
            obspower(tpwr2);
            obspwrf(tpwr2f);
            delay(GRADIENT_RES-MIN_EVENT_TIME);
            obl_shapedgradient(ss2_grad.name,ss2_grad.duration,ss2_grad.amp*sign1,0.0,0.0,NOWAIT);
            delay(ss2_grad.rfDelayFront);
            shaped_pulse(p2pat,p2,oph,rof1,rof2);
            delay(ss2_grad.rfDelayBack);
            delay(delay3);
            //delay(region_delay2);
        

        }


        
/*****************************************************180(2)chirp pulse*********************************************************************/

   /*     obsoffset(tof);
        obspower(tpwr902);
        pulse(p902,zero);
        delay(delay3);
	delay(delay1);*/

region_delay2=p3;
            //delay(region_delay2);
            delay(4e-6);         
	    delay(delay3);
            obspower(tpwr3);
            obspwrf(tpwr3f);
            delay(GRADIENT_RES-MIN_EVENT_TIME);
            obl_shapedgradient(ss3_grad.name,ss3_grad.duration,0.0,ss3_grad.amp*sign1,0.0,NOWAIT);
            delay(ss3_grad.rfDelayFront);
            shaped_pulse(p3pat,p3,oph,rof1,rof2);
            delay(ss3_grad.rfDelayBack);
            delay(delay3);
	    
/**************************************************************crush grad********************************************************************/
 if (strchr(crusher,'x')||strchr(crusher,'y'))
                    {
                        if (strchr(crusher,'x')||strchr(crusher,'y')) obl_shapedgradient(spoil_grad.name,spoil_grad.duration,0.5*gcrushx*spoil_grad.amp,1.2*gcrushy*spoil_grad.amp,0.5*gcrushy*spoil_grad.amp,WAIT);
                    }
 
        delay(delay1);

        delay(delay4);
/*************************************************90segment_selected_sinc pulse******************************************************/
        obsoffset(regiontof);
        obspower(p1_rf.powerCoarse);
        obspwrf(p1_rf.powerFine);
        delay(GRADIENT_RES-MIN_EVENT_TIME);
        delay(GRADIENT_RES-MIN_EVENT_TIME);

        obl_shapedgradient(regionr_grad.name,regionr_grad.duration,0,0,-regionr_grad.amp,WAIT);
        obl_shapedgradient(region_grad.name,region_grad.duration,0.0,0.0,region_grad.amp,NOWAIT);
        delay(region_grad.rfDelayFront);

        shaped_pulse(p1pat,region_grad.rfDuration,oph,rof1,rof2);
        delay(region_grad.rfDelayBack);


        obsoffset(tof);
        delay(delay3);



 /*********************************************************crush grad********************************************************************/
        //delay(region_delay3);
        if (strchr(crusher,'x')||strchr(crusher,'y'))
        {
            if (strchr(crusher,'x')||strchr(crusher,'y')) obl_shapedgradient(spoil_grad.name,spoil_grad.duration,gcrushx*spoil_grad.amp,gcrushy*spoil_grad.amp,gcrushy*spoil_grad.amp,WAIT);
        }
 

    /* Begin compressed multislice loop *********************/
        for(rrr=0;rrr<ns;rrr++)
        {
            if (EPI && !trtype) delay(delay1); //delay(tr_delay); /* EPI relaxation delay */

            ifzero(vprep); /* if preparation modules are to be applied */

            sp1on(); delay(4e-6); sp1off(); /* scope trigger */

            endif(vprep); /* end of preparation modules */

            for (j=0;j<ncseg;j++)
            { /* loop over compressed segments */
                if (CSEG && RSEG)
                { /* Reverse reference segment order */}
                else l=j;

                if (ALT)
                { /* Alternate read gradient on alternate segments */
                    if ((i+j)%2) sign = -1.0;
                    else sign = 1.0;
                }
                if (EPI)
                { /* EPI */
/*************************************************90slice_selected_sinc pulse******************************************************/
                    slice_thk=thk;
                    slicetof=tof+(region/slice_thk)*(regiontof-tof)-425.8*slice_thk*slice_grad.amp*(ns-1)/2+425.8*slice_thk*slice_grad.amp*rrr;
                    
                    //if (regiontof<tof) slicetof=tof+(region/slice_thk)*(regiontof-tof)+425.8*slice_thk*slice_grad.amp*(ns-1)/2-425.8*slice_thk*slice_grad.amp*rrr;
                    //else slicetof=tof+(region/slice_thk)*(regiontof-tof)-425.8*slice_thk*slice_grad.amp*(ns-1)/2+425.8*slice_thk*slice_grad.amp*rrr;
                    //if (regiontof<tof) slicetof=tof+(region/slice_thk)*(regiontof-tof)+425.8*slice_thk*slice_grad.amp*(ns+1)/2-425.8*slice_thk*slice_grad.amp*(rrr+1);
                    //else slicetf=tof+(region/slice_thk)*(regiontof-tof)-425.8*slice_thk*slice_grad.amp*(ns+1)/2+425.8*slice_thk*slice_grad.amp*(rrr+1);
                    obsoffset(slicetof);
                    obspower(p1_rf.powerCoarse);
                    obspwrf(p1_rf.powerFine);
                    delay(GRADIENT_RES-MIN_EVENT_TIME);
                }

                if (SE || DW || VS)
                { /* Spin echo EPI or EPI with diffusion or vascular suppresion */
                    obl_shapedgradient(slice_grad.name,slice_grad.duration,0.0,0.0,slice_grad.amp,NOWAIT);
                    delay(slice_grad.rfDelayFront);
                    shapedpulse(p1pat,slice_grad.rfDuration,oph,rof1,rof2);
                    delay(slice_grad.rfDelayBack);
                    obl_shapedgradient(slicer_grad.name,slicer_grad.duration,0,0,-slicer_grad.amp,WAIT);
                    delay(8e-6);
                    obsoffset(tof);

/**************************************************************crush grad********************************************************************/
                    if (strchr(crusher,'x')||strchr(crusher,'y'))
                    {
                        if (strchr(crusher,'x')||strchr(crusher,'y')) obl_shapedgradient(spoil_grad.name,spoil_grad.duration,0.5*gcrushx*spoil_grad.amp,1.2*gcrushy*spoil_grad.amp,0.5*gcrushy*spoil_grad.amp,WAIT);
                    }
                    delay(delay2);
                   // obl_shapedgradient(ror_grad.name,ror_grad.duration,-roramp*sign*sign_ror,0,0,NOWAIT);
                   // delay(ror_grad.duration-alfa);
                    startacq(alfa); /* prepare for acquisition */
                }

                if (EPI)
                { /* EPI */
                    if (TC) delay(tc_delay[i+l]); /* segment time correction */
            		//obl_shaped3gradient(epirR.name,epipename[i+l],epirS.name,epirTime,roamp*sign*sign_ror,peamp[i+l]*sign_per,epirS.amp,NOWAIT);
                    //obl_shaped3gradient(ropat,pepat,epirS.name,epirTime,gradro,gradpe,epirS.amp,NOWAIT);
		        me=i%4;
            switch(me){
                case 0:obl_shaped3gradient(ropat1,pepat1,epirS.name,epirTime,gradro1,gradpe1,epirS.amp,NOWAIT);break;
                case 1:obl_shaped3gradient(ropat2,pepat2,epirS.name,epirTime,gradro2,gradpe2,epirS.amp,NOWAIT);break;
                case 2:obl_shaped3gradient(ropat3,pepat3,epirS.name,epirTime,gradro3,gradpe3,epirS.amp,NOWAIT);break;
                case 3:obl_shaped3gradient(ropat4,pepat4,epirS.name,epirTime,gradro4,gradpe4,epirS.amp,NOWAIT);break;
            }

        
            		nwloop(ntraces,vtraces,vtrace_ctr);
            		//delay(delacq);
			printf("samplenp is %d\n",samplenp);
            		printf("delacq is %.2f us\n",delacq*1e6);

/**************************************************************Sampling********************************************************************/
            		for(k=0;k<samplenp;k++)
            		{
                		sample(dw);
                		delay(si[k]);
            			
			}
            		delay(delacq);
            		endnwloop(vtrace_ctr);
            		if (aqtm > at) sample(aqtm-at); /* finish sampling (if required) */
            		endacq();    /* finalise the acquisition */
            		delay(4e-6); /* for very large sw it seems there can be a slight timing error ?? */
            		if (TC) delay(tcmax-tc_delay[i+l]);  /* segment time correction */
            		if (CSEG) delay(segtr_delay); /* compressed segment delay */
                }
            } /* end of compressed EPI segment loop */

/**************************************************************Spoil********************************************************************/
            if (EPI && tspoil>0.0)
            { /* compressed EPI segment spoiler */
                //obl_shapedgradient(spoil_grad.name,spoil_grad.duration,0,0,0,WAIT);
                obl_shapedgradient(spoil_grad.name,spoil_grad.duration,spoil_grad.amp,spoil_grad.amp,spoil_grad.amp,WAIT);
            }
           delay(5000e-6);
        }
    } /* end of standard segment loop */

    /* Inter-image delay */
    sub(ntrt,ct,vtrimage);
    decr(vtrimage);
    ifzero(vtrimage);
    delay(trimage);
    endif(vtrimage);
}
