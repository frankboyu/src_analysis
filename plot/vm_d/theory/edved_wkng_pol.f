        program test
        common/par/pi,pm,dm,vmm
        common/input/ephin,sphin,bphin
        open(11,file='input.txt')
        read(11,*)ephin
        read(11,*)sphin
        read(11,*)bphin
        write(6,*)'ephin=',ephin,' sphin=',sphin,' bphin=',bphin
        pi = acos(-1.0)
        pm    = 0.938279
        kvm  = 3 ! 1-rho meson ! 3 phi meson  ! 4 J/Psi meson
        in   = 1 ! initialize
        call 
     & edved(in,kvm,ei,q2,q0,epsl,t,crs0,crs,tcrs0,tcrs,pd,thd,pvm,thvm, !initialization
     & crs0_m1,crs0_0,crs0_1,crs_m1,crs_0,crs_1,bphin,sphin)
       
 
        ei = 0.0 ! electron energy (for ed--> evd)


        x = 0.1 ! Bjorken x (for ed--> evd)

        q2 =  0.0  !(if q2=0 - real photoproduction)
        epsl = 0.0

        q0 = ephin                                    ! photon energy
        s = -q2 + 2.0*pm*q0 + pm**2
        w = sqrt(s)
        t_min = t_minimum(q2,s)
*        print *,t_min
        do it = 10,200,1                              ! defines t
        t = -float(it)/100.0
        if(t_min.lt.t)goto 1
*        print *,t
        in = 0
        call 
     & edved(in,kvm,ei,q2,q0,epsl,t,crs0,crs,tcrs0,tcrs,pd,thd,pvm,thvm,
     & crs0_m1,crs0_0,crs0_1,crs_m1,crs_0,crs_1,bphin,sphin)

        pdt  = pd*sin(thd)          !transverse momentum of deuteron
        pvmt = pvm*sin(thvm)        !transverse momentum of vmeson
        pdz  = pd*cos(thd)          !z momentum of deuteron
        pvmz = pvm*cos(thvm)        !z momentum of vmesons
        zmom = pdz + pvmz           !checking z momentum conservation

*******************************************************************
*     checking that individual polarization cross section give the
*     same unpolarized  and tensor polarized cross sections
******************************************************************        
        crs0_chk = (crs0_m1+crs0_0+crs0_1)/3.0
        crs_chk  = (crs_m1 + crs_0 + crs_1)/3.0
        
        tcrs0_chk = (crs0_m1 + crs0_1 - 2.0*crs0_0)/3.0
        tcrs_chk  = (crs_m1 + crs_1 - 2.0*crs_0)/3.0
       
        
        write(6,*)"       ",-t,crs0,crs,tcrs0,tcrs
*        write(6,*)"checkig",-t,crs0_chk,crs_chk,tcrs0_chk,tcrs_chk
*        write(12,12)abs(t),crs0,crs
 12     format(2x,f6.3,2(2x,e10.3))
 1      continue
        enddo
        end

 
        

       subroutine 
     & edved(in,ivm,ei,q2,q0,epsl,t,crs0,crs,tcrs0,tcrs,pd,thd,pvm,thvm,
     & crs0_m1,crs0_0,crs0_1,crs_m1,crs_0,crs_1,bphin,sphin)
**************************************************************************
*  in   - parameter for initialization (1)-initialize (0)-compute - input
*  ivm  - 1 - rho meson,  3- phi meson
*  ei   - initial energy of electron in GeV -  input (it is not used for the case of photoproduction)
*  q2   - (-) four-momentum square of virtual photon Gev^2 - input (=0 for the case of photoproduction)
*  q0   - virtual photon energy in Gev  - input  (real photon energy in the case of photoproduction)
*  epsl - polarization factor of virtual photon (0 for the case of photoproduction)
*  t    - transferred t of the reaction Gev^2(negative) - input         
*  crs0 - Born term(no FSI) of d\sigma/domega_edE'dt nanobarn/Gev^2 - output
*  crs  - full d\sigma/dOmega_e,dE',dt nanobarn/Gev^2 - output
*
* when q2 = 0, code calculates photoproduction cross section 
* crs0 - Born term(no FSI) of d\sigma/dt nanobarn/Gev^2 - output
* crs  - full d\sigma/dt nanobarn/Gev^2 - output
*  tcrs0 - Born tern for tensor polarized deuteron
*  tcrs - full cross section for tensor polarized deuteron
*   
*  pd   - momentum of final deuteron - GeV/c -output
*  thd  - polar angle of deuteron momentum vs q in radians - output
*  pvm  - momentum of the vector mesons in GeV/c - output
*  thvm - polar angle of vector mesons vs q in radians - output
*
* edved.f version 01-June-2000, FIU       
* modified only for phi meson production
* with new version of parameterization of gamma N --> phi N amplitude
* 16-Aug-2006, FIU  
* rho - meson production is added to this version of the code
* 9-March - 2010
* J/Psi - meson production is added 
*     5-April-2012
*     Scattering from tensor polarized target is added
*     10-March-2025
*     Output now contains cross section for individual polarization of
*     the deuteron
*       sigma0_m1, sigma0_0, sigma0_1  for -1, 0, 1 polarizations for PWIA
*       sigma_m1, sigma_0, sigma_1  for -1, 0, 1 polarizations for Full
*     09/06/25
*     
*        
**************************************************************************
	common/par/pi,pm,dm,vmm
        common/alpha/alpha
      	common/photon/q2c,q0c,qv
        common/bjor/x
	common/photon_virt/qq2
        common/idepk/idepk
        common/ctornot/ict0
        common/cross_gn/sigma_gn
        common/realpart_gn/alpha_g
        common/slope_gn/b_g
        common/case_gn/icase
        common/slope_vn/b_vn
        common/cross_vn/sigma_vn
        common/real_vn/al_vn
        common/gn_scale/gn_scale
        common/vm_case/kvm

        if(in.eq.1)then				
********************************************************************
*         Initialization
********************************************************************
           kvm = ivm



*************** P A R A M E T E R S **************************************

**************** rho - Meson *********************************************
           if(ivm.eq.1)then     ! rho meson parameters


           vmm = 0.77           ! mass of rho meson   
****************************************************************************
* Parameters of gamma + N --> rho + N amplitude         *
****************************************************************************
*   gn_scale - rescales the gamma N -->  rho N cross section
*   b_g     - slope factor for thr gamma - N --> rho N diff cross section
*   alpha_g - real part of the  gamma - N --> rho N  amplitude
****************************************************************************
           gn_scale =   1.0
           b_g      =   8.0
*           b_g     = 8.7 ! for egm=6 GeV

           alpha_g  =  -0.4 !-0.27
***************************************************
* parameterization of dsigma/dt at t=0 given in the 
* "function dsdt_gn"   as sigma_gn_0 - 
*this may be updated for JLab energies.
***************************************************

 
****************************************************************************
* Parameters of  rho + N --> rho' + N' amplitude
****************************************************************************
*  sigma_vn  -  rho-N cross section in milibarns
*  b_vn      - slope factor of the amplitude
*  al_vn     - real part of the amplitude
****************************************************************************
        sigma_vn = 28
*        sigma_vn = 32 ! for egm=6 GeV
        b_vn     = 8.0
*        b_vn     = 8.7 ! for egm=6 GeV
        al_vn    = -0.4


**************** Phi Meson *****************************
          elseif(ivm.eq.3)then ! phi meson parameters


           vmm   = 1.0197     ! phi meson mass
********************************************************************
* Parameters of gamma + N --> phi + N amplitude         *
*********************************************************
*   gn_scale - rescales the gamma N -->  phi N cross section
*   icase =1,2,3,4  
*   (1 and 2) are the parameterizations from Eq. 3.85a and 3.85b 
*      from BYSP
*   (3) uses the parameterization of rhe curve of dsigma/dt (t=tmin)
*      given in Fig. 3 of T. Mibe at al Phys.Rev. Lett. 95 182001,2005
*      Note that for this case the  t dependence is parameterzied as
*      exp(b_g/2 t) where slope factor b_g is used as a free 
*      parameter. I did not use  3.38  quoted in T. Mibe et al  PRL, since
*      it is defined for restricted range of t and photon energies.
*   (4) is the hybrid of (2) and (3). At t_min the ds/dt is used as in icase=3
*      while the t dependence is parameterizied according to icase=2.
*      This allows to check the importance of t^2 term in the exponent 
*      of the amplitude  
*
*   b_g     - slope factor, relevant only for icase=3
*   alpha_g - real part of the amplitude - relevant for all icases

********************************************************************
        gn_scale = 1.0
        icase    =   1
        b_g      =   5.0 !relevant only for icase=3,4,5
        alpha_g  =  -0.5
********************************************************************
*         Parameters of cross section, slope  factor and real part
********************************************************************
* Parameters of phi + N --> phi' + N' amplitude
*
*  sigma_vn  -  phi-N cross section in milibarns
*  b_vn      - slope factor of the amplitude
*  al_vn     - real part of the amplitude
********************************************************************
        sigma_vn = sphin
        b_vn     = bphin
        al_vn    = -0.5
        write(6,*)'phi-N parameters: sigma_vn=',sigma_vn,' b_vn=',b_vn
********************************************************************

************* J/PSI ************************************************
        elseif(ivm.eq.4)then

           vmm = 3.096916   !  mass of J/Psi meson

****************************************************************************
* Parameters of gamma + N --> J/Psi + N amplitude         *
****************************************************************************
*   gn_scale - rescales the gamma N -->  J/Psi N cross section
*   b_g     - slope factor for thr gamma - N --> J/Psi N diff cross section
*   alpha_g - real part of the  gamma - N --> J/Psi N  amplitude
****************************************************************************
           gn_scale =   1.0
           sigma_gn = 1.01E-3 !microbarn/GeV^2
           b_g      =   1.25
           alpha_g  =  -0.5
  
****************************************************************************
* Parameters of  J/Psi + N --> J/Psi' + N' amplitude
****************************************************************************
*  sigma_vn  -  J/Psi-N cross section in milibarns
*  b_vn      - slope factor of the amplitude
*  al_vn     - real part of the amplitude
****************************************************************************
        sigma_vn = 4.0  
        b_vn     = 2.0
        al_vn    = -0.5

****************************************************************************
        endif

        ict = 0
        ins = 1
        call form_factors(ict,ins,qt,qz,q2,fc,fq,tc,tq) 
        pi    = acos(-1.0)
        alpha = 1.0/137.0
        pm    = 0.938279
        dm    = 1.875628
        return

        else 

****************************************
*        Calculation
***************************************
        ict = 0
********************************************************************
*         Input parameters
********************************************************************
        q0c = q0
        q2c = q2
        qq2 = q2
        x = q2/2.0/pm/q0
        qv = sqrt(q2 + q0**2)
        er = ei - q0 
*        s  = pm**2 + 2.0*pm*q0 - q2
        
*******************************************************
*                         Zroyacum
*******************************************************
 			         sigma0 = 0.0
			         sigma = 0.0

*******************************************************
* Kinematics of coherent \rho electroproduction
*      q + p_d --> p_\rho + p_d'
*******************************************************
	q_d_min = (vmm**2+Q2)/2.0/qv
        t0      = t
*	t_min   = q_d_min**2
*        t_min = t_minimum(q2,s)
*        if(abs(t0).lt.t_min)t = -t_min

 	q_d     = sqrt(-t + t**2/4.0/dm**2)
	e_d_pr  = sqrt(dm**2 + q_d**2) 
	q_d_z   = (vmm**2+q2-t-2.0*q0*(dm-e_d_pr))/2.0/qv

*        if(abs(t0).lt.t_min)then
*        q_d_z = (vmm**2+Q2)/2.0/qv
*        q_d   = q_d_z
*      	 e_d_pr  = sqrt(dm**2 + q_d**2) 
*        endif

	th_d_pr = 0.0
	if(q_d.ne.0.0)then
	arg =   (q_d_z / q_d)
	if(arg**2.gt.1.0)arg=0.0 
	th_d_pr =  acos(arg)
	endif
	qz = -(q_d_z - (e_d_pr - dm))


        arg2 = q_d**2 - q_d_z**2
        if(arg2.le.0.0)arg2=0.0
        q_d_t = sqrt(arg2)


*ds/domega/de' dt
********************************************************
 	w2    = -q2 + 2.0*q0*pm + pm**2
	eff_k = (w2-pm**2)/2.0/pm
	g_v1  = alpha/4.0/pi**2  
        g_v = 1.0
        
        if(q2.gt.0.0)
     &	g_v   = g_v1 * eff_k/q2 * er/ei * 2.0/(1.0-epsl) 
        R     = 0.4
	rlt   = R * (q2/vmm**2)
        fct   = 1.0/(1.0+q2/vmm**2)**2

*        write(6,*)g_v,rlt

        ict0 = ict 
        ita  = 1
        imu = 1
	call sigma_0(ita,imu,sh0_1,sh1_1,sh2_1,sh_1,t,q_d_t,qz,
     &                                                    q_d_z)
*        print *,ita,imu,sh0_1,sh1_1,sh2_1,sh_1,t,s,q_d_t,qz,q_d_z
        imu = 0
	call sigma_0(ita,imu,sh0_0,sh1_0,sh2_0,sh_0,t,q_d_t,qz,
     &                                                    q_d_z)
        imu = -1
	call sigma_0(ita,imu,sh0_m1,sh1_m1,sh2_m1,sh_m1,t,q_d_t,qz,
     &                                                    q_d_z)

        sh0 = (sh0_1 + sh0_0 + sh0_m1)/3.0 
        sh1 = (sh1_1 + sh1_0 + sh1_m1)/3.0 
        sh2 = (sh2_1 + sh2_0 + sh2_m1)/3.0 
        sh  = (sh_1 + sh_0 + sh_m1)/3.0 

*        print *, "aha",sh0,sh1,sh2

*****************************************************
*     Calculation from tensor polarized deuteron
******************************************************        
        th0 = (sh0_1 + sh0_m1- 2.0*sh0_0 )/3.0 
        th  = (sh_1  + sh_m1 - 2.0*sh_0 )/3.0 
*************************************************************
*       Unpolarized case
*************************************************************
	s_t0 = sh0 * fct 
	s_t  = sh  * fct 

	s_l0 = s_t0 * rlt
	s_l  = s_t  * rlt

	sigma0 = g_v * (s_t0 + epsl*s_l0)
	sigma  = g_v * (s_t  + epsl*s_l )

        crs0 = sigma0
        crs  = sigma
***********************************************************
*     06/09/25
*     Calculation of the cross section for given polarization
*     of the deuteron
***********************************************************
        s_t0_m1 = sh0_m1  * fct
        s_t0_0  = sh0_0   * fct 
        s_t0_1  = sh0_1   * fct
        
	s_l0_m1 = s_t0_m1 * rlt
	s_l0_0  = s_t0_0  * rlt
	s_l0_1  = s_t0_1  * rlt
        
	sigma0_m1 = g_v * (s_t0_m1 + epsl*s_l0_m1)
	sigma0_0  = g_v * (s_t0_0  + epsl*s_l0_0)
	sigma0_1  = g_v * (s_t0_1  + epsl*s_l0_1)

        
        s_t_m1 = sh_m1 * fct
        s_t_0  = sh_0  * fct 
        s_t_1  = sh_1  * fct
        
	s_l_m1 = s_t_m1 * rlt
	s_l_0  = s_t_0  * rlt
	s_l_1  = s_t_1  * rlt

 	sigma_m1 = g_v * (s_t_m1 + epsl*s_l_m1)
	sigma_00  = g_v * (s_t_0  + epsl*s_l_0)
	sigma_1  = g_v * (s_t_1  + epsl*s_l_1)
              

        crs0_m1 = sigma0_m1
        crs0_0  = sigma0_0
        crs0_1  = sigma0_1
        
        crs_m1 = sigma_m1
        crs_0  = sigma_00
        crs_1  = sigma_1

        
        
*************************************************************
*       Tensor Polarized Case
*************************************************************
	th_t0 = th0 * fct 
	th_t  = th  * fct 

	th_l0 = th_t0 * rlt
	th_l  = th_t  * rlt

	tsig0 = g_v * (th_t0 + epsl*th_l0)
	tsig  = g_v * (th_t  + epsl*th_l )

        tcrs0 = tsig0
        tcrs  = tsig

**************************************************
* scattered deuteron kinematics
        pd  = q_d
        thd = th_d_pr

**************************************************
*      Scattered vector mesons kinematics
**************************************************	
	p_rho   = sqrt(qv**2 - 2.0*qv*q_d_z + q_d**2)
	p_rho_z = qv - q_d_z 
	th_rho = 0.0
	if(p_rho.ne.0.0)then
	arg = p_rho_z / p_rho 
	if(arg.gt.1.0) arg = 1.0
	if(arg.lt.-1.0)arg = -1.0
	th_rho =  acos(arg)
	endif
        pvm  = p_rho
        thvm = th_rho

        return
        endif        
 101    end


       function t_minimum(q2,s)
        implicit none
        common/par/pi,pm,dm,vmm
        real pi,pm,dm,vmm
        real t_minimum
        real q2,s
        real en_cm,pn_cm
        real ev_cm,pv_cm
        en_cm = (s + q2 + pm**2)/(2.0*sqrt(s))
        pn_cm = sqrt(en_cm**2 - pm**2)
        ev_cm = (s - pm**2 + vmm**2)/(2.0*sqrt(s))
        pv_cm = sqrt(ev_cm**2 - vmm**2)
*        print *,"*",q2,s
        t_minimum = -s - q2 + 2.0*(ev_cm*en_cm + pv_cm*pn_cm) + pm**2
        return
        end


 	subroutine sigma_0(ita,imu,sh0,sh1,sh2,sh,t,qt,qz,qdz)
*********************************************************
*       VMD model for \sigma_T (electroproduction)
*********************************************************
 	common/par/pi,pm,dm,vmm
      	common/photon/q2,q0,qv
        common/idepk/itb
        common/kvm/kvm
*********************************************************

*********************************************************
        fct1 = 1.0/(16.0*pi)
	fct2 =  0.389385 *1000.*1000.
*********************************************************
        qx  = qt
        qy  = 0.0
        p_g = qv  !(s-pm**2)/(2.0*pm)

        sh0 = 0.0
        sh0 = fct1*FaFa(ita,imu,p_g,t,qx,qy,qz)

        qpz = -qdz/2.0 - t/8.0/pm + (q2+vmm**2)/2.0/qv

        sh1 = 0.0
        sh1 = fct1*two_FaFb(ita,imu,p_g,t,qx,qy,qz,qpz)

        qppz = -qdz/2.0 - t/8.0/pm + (q2+vmm**2)/2.0/qv
        
        sh2 = 0.0        
        sh2 = fct1*FbFb(ita,imu,p_g,t,qx,qy,qz,qpz,qppz)

	sh0 = sh0 * fct2 
	sh1 = sh1 * fct2 
	sh2 = sh2 * fct2 
*        print *,sh0,sh1,sh2
 	sh = sh0 + sh1 + sh2

 	return
	end




        function FaFa(ita,imu,p_g,t,qx,qy,qz)
        common/par/pi,pm,dm,vmm
      	common/photon/q2,q0,qv

        qt   = sqrt(qx**2 + qy**2) 
        qq   = sqrt(qx**2 + qy**2 + qz**2)
        s    = pm**2 + 2.0*pm*q0 - q2
*        t     = -qq**2
        a1   = f_gN(s,t)**2
        a2   = (1.0+ al_gN(s)**2)*a1
        qxk  = qx/2.0
        qyk  = qy/2.0
        qzk  = qz/2.0
        ms   = 1 
        den  = ro(ms,1,1,ita,imu,qxk,qyk,qzk,qxk,qyk,qzk)
        FaFa = 4.0*a2*den 
*        print *,"fafa",fafa,"s=",s,t,a1
        return
        end

        function two_FaFb(ita,imu,p_g,t,qx,qy,qz,qpz)
        external under_ab,phip_a,phip_b
        common/for_ab/iita,iimu,pp_g,tt,qqx,qqy,qqz,qqpz
        common/par/pi,pm,dm,vmm
        iita = ita
        iimu = imu
        pp_g = p_g
        tt   = t
        qqx  = qx
        qqy  = qy
        qqz  = qz
        qqpz = qpz
        qp_a = 0.0
        qp_b = 1.6
        eps = 0.001
        call gadap2(qp_a,qp_b,phip_a,phip_b,under_ab,eps,sum)
        two_FaFb = sum/((2.0*pi)**2)
        return
        end 

        function under_ab(qpt,phip)
        common/for_ab/ita,imu,p_g,t,qx,qy,qz,qpz
        common/par/pi,pm,dm,vmm
      	common/photon/q2,q0,qv    
        common/ctornot/ict

        s = pm**2 + 2.0*pm*q0 - q2
        qpx   = qpt*cos(phip)
        qpy   = qpt*sin(phip)

	q_t     = sqrt(qx**2 + qy**2)
	qqp_t   = qx*qpx  + qy*qpy  

	q       = sqrt(qx**2   + qy**2  + qz**2)
	qp      = sqrt(qpx**2  + qpy**2 + qpz**2)
	qqp     = qx*qpx  + qy*qpy  + qz*qpz
*        t = - q**2

        arg_m = q**2/4.0 - qqp + qp**2        
        if(arg_m.lt.0.0)argm=0.0
        t_m   = -arg_m
	q_m  = sqrt(arg_m)


        arg_p = q**2/4.0 + qqp + qp**2
        if(arg_p.lt.0.0)arg_p = 0.0 
        t_p   = -arg_p
	q_p  = sqrt(arg_p)
        

        t_1   = f_gn(s,t)*f_gn(s,t_m)*f_vn(s,t_p)
        a_1   = (1 + al_gn(s)*al_gn(s))

        am_re = -a_1*t_1
        am_im = al_vn(s)*am_re
        
        qkx = qx/2.0
        qky = qy/2.0
        qkz = qz/2.0
                   ms = 1
        ro_re = ro(ms,1,1,ita,imu,qkx,qky,qkz,qpx,qpy,qpz)
        ro_im = ro(ms,2,1,ita,imu,qkx,qky,qkz,qpx,qpy,qpz)

        under_ab_1 = 2.0*(am_re*ro_re - am_im*ro_im)


        am_re_d = al_vn(s)*a_1*t_1 
        am_im_d = - a_1 * t_1

        ro_re_d = ro(ms,1,2,ita,imu,qkx,qky,qkz,qpx,qpy,qpz)
        ro_im_d = ro(ms,2,2,ita,imu,qkx,qky,qkz,qpx,qpy,qpz)

        term = - 4.0/sqrt(2.0*pi)

        under_ab_2 = term*(am_re_d*ro_re_d - am_im_d*ro_im_d)

        under_ab =  (under_ab_1 + under_ab_2)*qpt

        return
        end


        function phip_a(x)
        phip_a  = 0.0
        return
        end

        function phip_b(x) 
        phip_b  = 2.0*acos(-1.)
        return
        end



        function FbFb(ita,imu,p_g,t,qx,qy,qz,qpz,qppz)
        external under_bb,phip_a,phip_b
        common/for_bb/iita,iimu,pp_g,qqx,qqy,qqz,qqpz,qqppz
        common/par/pi,pm,dm,vmm
        iita  = ita
        iimu  = imu
        pp_g  = p_g
        qqx   = qx
        qqy   = qy
        qqz   = qz
        qqpz  = qpz
        qqppz = qppz        
        qp_a = 0.0
        qp_b = 1.6
        eps = 0.001
        call gadap2(0.0,1.6,phip_a,phip_b,under_bb,eps,sum)
        FbFb = sum/((2.0*pi)**2)/((2.0*pi)**2) *100.0
        return
        end 

        function under_bb(qp,phip)
        external un_un_bb,phipp_a, phipp_b
        common/for_under_bb/qqpx,qqpy
        qqpx   = qp*cos(phip)
        qqpy   = qp*sin(phip)
        epsp = 0.001
        call gadaps2(0.0,1.6,phipp_a,phipp_b,un_un_bb,epsp,sm2)
        under_bb = sm2 *qp

        return
        end

        function phipp_a(x)
        phipp_a = 0.0
        return
        end

        function phipp_b(x)
        phipp_b = 2.0*acos(-1.0)
        return
        end


        function un_un_bb(qppt,phipp)
        common/for_bb/ita,imu,p_g,qx,qy,qz,qpz,qppz
        common/for_under_bb/qpx,qpy
        common/par/pi,pm,dm,vmm
      	common/photon/q2,q0,qv
        common/ctornot/ict
        
        s       = pm**2 + 2.0*pm*q0 - q2
 	q_t     = sqrt(qx**2 + qy**2)
 	q       = sqrt(qx**2 + qy**2 + qz**2)


	qp_t    = sqrt(qpx**2  + qpy**2)
	qqp_t   = qx*qpx  + qy*qpy
	qp      = sqrt(qpx**2  + qpy**2 + qpz**2)
	qqp     = qx*qpx  + qy*qpy + qz*qpz


	q_m  = sqrt(q**2/4.0 - qqp + qp**2)
        t_m  = -q_m**2

	q_p  = sqrt(q**2/4.0 + qqp + qp**2)
        t_p  = -q_p**2


        qppx = qppt*cos(phipp)
        qppy = qppt*sin(phipp)

	qqpp_t   = qx*qppx  + qy*qppy 
	qpp    = sqrt(qppx**2  + qppy**2 + qppz**2) 
	qqpp   = qx*qppx  + qy*qppy + qz*qppz 

	qq_m  = sqrt(q**2/4.0 - qqpp + qpp**2)
        tt_m  = -qq_m**2
	qq_p  = sqrt(q**2/4.0 + qqpp + qpp**2)
        tt_p  = -qq_p**2

*        print *,"t_m",t_m,tt_m
        t_1   = f_gn(s,t_m)*f_vn(s,t_p)
        t_2   = f_gn(s,tt_m)*f_vn(s,tt_p)

        a_gn = (1 + al_gn(s)*al_gn(s))
        a_vn = (1 + al_vn(s)*al_gn(s))
                    ms = 1
        ro_re    = ro(ms,1,1,ita,imu,qpx,qpy,qpz,qppx,qppy,qppz)
        ro_im_d  = ro(ms,2,2,ita,imu,qpx,qpy,qpz,qppx,qppy,qppz)
        ro_re_dd = ro(ms,1,3,ita,imu,qpx,qpy,qpz,qppx,qppy,qppz)

        un_un_bb_1 =  1.0/4.0* a_gn*a_vn*t_1*t_2*ro_re
        un_un_bb_2 =  1.0/sqrt(2.0*pi)*a_gn*a_vn*t_1*t_2*ro_im_d
        un_un_bb_3 =  1.0/(2.0*pi)*a_gn*a_vn*t_1*t_2*ro_re_dd

        un_un_bb = (un_un_bb_1  + un_un_bb_2 +  un_un_bb_3)*qppt /100.00

        return
        end 


        function f_gn(s,t)
**************************************************
* Imagniary part of  gamma N --> vm  N amplitude
**************************************************
      common/case_gn/icase
      common/par/pi,pm,dm,vmm
      common/gn_scale/gn_scale
      common/vm_case/kvm
      scale = 1.0
      if(kvm.eq.1)then
      call f_gN_vmN(s,t,f_re,f_im,icase,kvm)

      elseif(kvm.eq.3)then 
       if(icase.eq.4)then
       eg = (s-pm**2)/2.0/pm
       tmin = -(vmm**2/(2.0*eg))**2.
       call  f_gN_vmN(s,tmin,f_re,f_im3,3,kvm)
       call  f_gN_vmN(s,tmin,f_re,f_im2,2,kvm)
       scale = f_im3/f_im2  
       endif
       call f_gN_vmN(s,t,f_re,f_im,icase,kvm)
       elseif(kvm.eq.4)then
*          print *,"kvm=",kvm,s,t,icase
       call f_gN_vmN(s,t,f_re,f_im,icase,kvm)
*       print *,"f_im=",f_im
      endif

      f_gn  = f_im*scale*gn_scale
*      print *,"scale",eg,scale
      return
      end



      subroutine f_gN_vmN(s,t,f_re,f_im,icase,kvm)
********************************************************************
*  Real and Imaginary parts of the gamma N--> vm N amplitude
*  are in Gev^{-2}
********************************************************************
      common/par/pi,pm,dm,vmm
      al = alg(s,kvm)
      bngev = 389.379
      term0 = dsdt_gn(s,t,icase,kvm)/bngev*(16.0*pi)/(1.0+al**2)
*      print *,"term0=",term0,al,pi,s,t,icase,kvm
      term  = sqrt(term0)
      f_re = term*al
      f_im = term
      return
      end

 
      function dsdt_gn(s,t,icase,kvm)
*******************************************************************
*  This is a parameterization for the differential
*  cross section of gamma+p-> vm + p
* dsigma/dt is in microbarn/GeV^2
* s = m^2 + 2.0*pm*e_gamma  - Q^2 in the lab frame
* 
*
********************************************************************
      common/par/pi,pm,dm,vmm
      common/cross_gn/sigma_gn
      common/photon_virt/q2
      dsdt_gn = 0.0
      eg = (s - pm**2)/(2.0*pm)
*      tmin = - (vmm**2/(2.0*eg))**2
      tmin =  t_minimum(q2,s)
      tpr = t-tmin
      dsdt = 0.0              
      if(tpr.gt.0.0)return

      if(kvm.eq.1)then 
***************for Rho Meson ************************************
* Paremeterizations are from BSYP
* 9-March-2010
*******************************************************************
      sigma_gn_0 = 563.87 - 259.58*eg + 57.714*eg**2 
     &            -5.7068*eg**3 + 0.20808*eg**4
      
      if(eg.gt.5.9)sigma_gn_0 = 110.0
      if(eg.gt.10.0)sigma_gn_0 = 90.0
      dsdt = sigma_gn_0* exp(b_gn(s,kvm)*t)
********************************************************************

      elseif(kvm.eq.3)then 
************** for Phi Meson ************************************
* Paremeterizations are from BSYP, page 340 (icase=1,2)
*                   and from Mibe, PRL 2005,3 Fig (icase=3)
*
*
* correction
* icase=3 case  when |t| is less than |tmin| cross section is set 0 
* Misak Sargsian
* 13-Aug-06
* Miami
********************************************************************
       if(icase.eq.1)then
       dsdt = 2.59*exp(5.9*t+ 1.4*t**2)
       elseif(icase.eq.2.or.icase.eq.4)then
       alpha = 1.14 + 0.27*t
       dsdt  = s**(2.*(alpha-1.))*1.34*exp(4.8*t + 1.7*t**2)
       elseif(icase.eq.3)then
                    tpr = t-tmin
                    dsdt = 0.0              
      !if(tpr.le.0.0)
       dsdt = dsdt_tmin(eg)*exp(b_gn(s,kvm)*tpr)
       elseif(icase.eq.5)then
       dsdt = sigma_gn*exp(b_gn(s,kvm)*t)
       endif
************************************************************************

       elseif (kvm.eq.4)then
       dsdt = sigma_gn*exp(b_gn(s,kvm)*t)
*       dsdt = 1.01E-3*exp(b_gn(s,kvm)*t)

*       print *,"dsdt",dsdt,b_gn(s,kvm)

      endif       
      if(dsdt.le.0.0)dsdt=0.0

 1     dsdt_gn = dsdt
      return
      end



      function dsdt_tmin(eg)
****************************************************************
* This fit is for the eg -- 1.6-6.6 GeV
* This fits Fig.3 of Mibe et al PRL (nucl-ex/0506015, 2005)
* scaled by 0.85
****************************************************************
      dsdt_tmin =  -6.9937  +8.2732*eg - 3.2595*eg**2
     &          +0.67068*eg**3 -0.69670E-01*eg**4 +0.28879E-02*eg**5   
      dsdt_tmin = dsdt_tmin
      return
      end

  

        function al_gn(s)
        common/vm_case/kvm
        al_gn = alg(s,kvm)
        return
        end
        function alg(s,kvm)
        common/realpart_gn/alpha_g
        alg = alpha_g
        return
        end



        function b_gn(s,kvm)
        common/slope_gn/b_g
        b_gn = b_g
        return
        end


        function f_vn(s,t)
***************************************************
* Function parameterizes vector-meson nucleon 
* scattering amplitude
***************************************************
	common/par/pi,pm,dm,vmm
        common/slope_vn/b_vn
        common/cross_vn/sigma_vn_0
*        sigma_vn_0 = 15.0
*        b_vn = 5.6
*        f_vn = 0.0
*        ep  = (s - pm**2 - vmm**2)/(2.0*pm)
*        if(ep.le.vmm)return
*        p   = sqrt(ep**2 - vmm**2)
        sigma_vn = sigma_vn_0 / 0.389385 
        f_vn   = sigma_vn * exp(b_vn/2.0*t)
        return
        end

        function al_vn(s)
        common/real_vn/an
*         an = -0.5
        al_vn = an
        return
        end

      



        function ro(ms,ir,iro,ita,imu,q1x,q1y,q1z,q2x,q2y,q2z)
	common/photon_virt/qq2
*******************************************************************************
*        ir  = 1 -ReR,  2 -ImR
*        iro = 1 -R,    2 - Delta R  3 - Delta Delta R
*        ita = 1 - polarizatio v.s. of virtual photon or final vector meson k
*            = 2 - vs of transfered momentun-q , 3 - vs [kxq] n
*        imu = 0,1-1 projection number
*******************************************************************************
        q1   = sqrt(q1x**2 + q1y**2 + q1z**2)
        q1t  = sqrt(q1x**2 + q1y**2)
        q2   = sqrt(q2x**2 + q2y**2 + q2z**2)
        q2t  = sqrt(q2x**2 + q2y**2)
        q1q2 = q1x*q2x + q1y*q2y + q1z*q2z

*	if(ms.eq.1)then
	ict1=0
	ict2=0

        if(iro.eq.1)then      ! R function  
        call form_factors(ict1,0,q1t,q1z,qq2,fc_1,fq_1,tc,tq) 
        call form_factors(ict2,0,q2t,q2z,qq2,fc_2,fq_2,tc,tq) 
        fc1  = fc_1
        fc2  = fc_2
        fq1  = fq_1
        fq2  = fq_2
        elseif(iro.eq.2)then  ! Delta R function
        call form_factors(ict1,0,q1t,q1z,qq2,fc_1,fq_1,tc,tq) 
        call form_factors(ict2,0,q2t,q2z,qq2,fc_2,fq_2,tc_2,tq_2) 
        fc1  = fc_1
        fc2  = tc_2
        fq1  = fq_1
        fq2  = tq_2
        elseif(iro.eq.3)then  ! Double Delta R function
        call form_factors(ict1,0,q1t,q1z,qq2,fc_1,fq_1,tc_1,tq_1) 
        call form_factors(ict2,0,q2t,q2z,qq2,fc_2,fq_2,tc_2,tq_2)  
        fc1  = tc_1
        fc2  = tc_2
        fq1  = tq_1
        fq2  = tq_2
        endif


        r1 = fc1*fc2

        if(ita.eq.1)then      ! tarberak I
        q1_1   = q1z
        q2_1   = q2z
        q1_2   = q1x**2 + q1y**2   
        q2_2   = q2x**2 + q2y**2
        q_12   = q1x*q2x + q1y*q2y  
        q_ir   = q1y*q2x - q1x*q2y
        elseif(ita.eq.2)then  ! tarberak II
        q1_1   = q1x
        q2_1   = q2x
        q1_2   = q1z**2 + q1y**2   
        q2_2   = q2z**2 + q2y**2
        q_12   = q1z*q2z + q1y*q2y
	q_ir   = q1y*q2z - q1z*q2y
        elseif(ita.eq.3)then  ! tarberak III
        q1_1   = q1y
        q2_1   = q2y
        q1_2   = q1x**2  + q1z**2   
        q2_2   = q2x**2  + q2z**2
        q_12   = q1x*q2x + q1z*q2z
	q_ir   = q1z*q2x - q1x*q2z
        endif

	if(ir.eq.1)then                    !real part 
        if(imu.eq.0)then
        r2m = 0.0
        if(q2.ne.0.0)r2m = 3.0*q2_1**2/q2**2-1.0
        r2 = r2m *fc1*fq2/sqrt(2.0)

        r3m = 0.0
        if(q1.ne.0.0)r3m = 3.0*q1_1**2/q1**2-1.0
        r3 = r3m *fc2*fq1/sqrt(2.0)

        r4m = 0.0
        if(q1.ne.0.0.and.q2.ne.0.0)then
        r4m = 9.0/2.0*q1_1*q2_1*q1q2/q1**2/q2**2
        r4m = r4m - 3.0/2.0*q1_1**2/q1**2 - 3.0/2.0*q2_1**2/q2**2
        r4m = r4m + 1.0/2.0
        endif
        r4  = r4m * fq1 * fq2

        else  ! for helicities +1 and -1
 
        r2m = 0.0
        if(q2.ne.0.0)r2m = 3.0*q2_2/q2**2/2.0 - 1.0
        r2 = r2m *fc1*fq2/sqrt(2.0)

        r3m = 0.0
        if(q1.ne.0.0)r3m = 3.0*q1_2/q1**2/2.0 - 1.0
        r3 = r3m *fc2*fq1/sqrt(2.0)

        r4m = 0.0
        if(q1.ne.0.0.and.q2.ne.0.0)then
        r4m = 9.0/2.0*q_12*q1q2/q1**2/q2**2/2.0
        r4m = r4m - 3.0/2.0*q1_2/q1**2/2.0 
        r4m = r4m - 3.0/2.0*q2_2/q2**2/2.0 
        r4m = r4m + 1.0/2.0
        endif
        r4  = r4m * fq1 * fq2
        endif

	elseif(ir.eq.2)then         !imaginary part
	r1 = 0.0

        if(imu.eq.0)then
        r2 = 0.0
        r3 = 0.0
        r4 = 0.0

        elseif(imu.eq.1)then  ! for helicities +1 and -1
 
        r2 = 0.0
        r3 = 0.0
        r4m = 0.0
        if(q1.ne.0.0.and.q2.ne.0.0)then
        r4m = 9.0/2.0*q_ir*q1q2/q1**2/q2**2/2.0
        endif
        r4  = r4m * fq1 * fq2
        elseif(imu.eq.-1)then  ! for helicities +1 and -1
        r2 = 0.0
        r3 = 0.0
        r4m = 0.0
        if(q1.ne.0.0.and.q2.ne.0.0)then
        r4m = 9.0/2.0*q_ir*q1q2/q1**2/q2**2/2.0
        endif
        r4  = - r4m * fq1 * fq2
        endif
	endif

        ro  = r1 + r2 + r3 + r4
        return
        end



      subroutine form_factors(ict,ins,qt,qz,q2,fc,fq,tc,tq) 
      common/formfactors/f_c(400),f_q(400),t_c(40,40),t_q(40,40)
      common/ctornot/ict0
      if(ins.eq.1)then
      open(unit=11,status='old',file='fc_fq.data')
      read(11,10)(f_c(k),f_q(k),k=1,400)
      close(11)
      open(unit=12,status='old',file='tc_tq.data')
*      write(12,10)((TC(kt,kz),TQ(kt,kz),kt=1,400),kz=1,400)
      read(12,10)((t_c(kt,kz),t_q(kt,kz),kt=1,40),kz=1,40)
      close(12)

*      open(unit=13,status='old',file='fc_fq_ct_mn.data')
*      read(13,10)(((f_c_ctmn(kq2,kt,kq),f_q_ctmn(kq2,kt,kq),
*     &               kq2=1,10),kt=1,11),kq=1,400)
*      close(13)

*      open(unit=14,status='old',file='fc_fq_ct_mx.data')
*      read(14,10)(((f_c_ctmx(kq2,kt,kq),f_q_ctmx(kq2,kt,kq),
*     &               kq2=1,10),kt=1,11),kq=1,400)
*      close(14)

10    format(8e10.3)  
      return
      endif
      q = sqrt(qt**2 + qz**2)
      i = q/0.005 + 1
      if(i.gt.400)i=400

      it = qt/0.04 + 1
      if(it.gt.40)it=40
      
      fct = 1.0
      if(qz.lt.0.0)fct = -1.0
      iz = abs(qz)/0.04 + 1
      if(iz.gt.40)iz=40

      iq2 = q2 + 0.1
      if(iq2.lt.1)iq2=1
      if(iq2.gt.10)iq2=10

      itt = qt**2*10 + 1.1
      if(itt.gt.11)itt=11
      itt = 1
 
*      if(ict.eq.0)then
*      fc = ffc(q)
*      fq = ffq(q)
      fc = f_c(i)
      fq = f_q(i)
*      fc = f_c(i)
*      fq = f_q(i)
      tc = fct*t_c(it,iz)
      tq = fct*t_q(it,iz)
      return
      end


      function ffc(q)
      if(q.le.0.25)then
      ffc =   0.99940     + 0.40200*q    - 93.710*q**2 + 690.44*q**3           
     &      -2128.3*q**4  + 2489.8*q**5 
      elseif(q.gt.0.25.and.q.le.0.6)then
      ffc =    1.1319    - 7.3529  *q    + 18.382*q**2 - 21.261*q**3           
     &       + 9.5674*q**4 
      elseif(q.gt.0.6)then
      ffc =    0.031985 - 0.23133 *q + 0.36184*q**2 - 0.20350*q**3           
     &       + 0.038222*q**4 
      endif
      if(ffc.gt.1)ffc=1.0
      return
      end

      function ffq(q)
      if(q.eq.0.0)then
      ffq = 0
      elseif(q.gt.0.0.and.q.le.0.25)then
      ffq =    0.14808E-03 - 0.40721E-01*q + 8.3114*q**2 -36.508*q**3  
     &       -91.903*q**4 + 866.69*q**5   - 1455.0*q**6
      elseif(q.gt.0.25.and.q.le.1.0)then
      ffq =    -0.44836E-01  + 1.3419*q - 5.4894*q**2 +10.115*q**3     
     &         -10.003*q**4  + 5.2052*q**5   - 1.1242*q**6
       elseif(q.gt.1.0)then
      ffq =     0.22920  - 0.65778*q + 0.73903*q**2 -0.40953*q**3     
     &        + 0.11248*q**4  - 0.12289E-01 *q**5   
      endif
      return
      end


                             
*********************************                                       
*   Subroutine for integration                                          
********************************                                        
      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAP1/ NUM,IFU                                            
      EXTERNAL F                                                        
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
    1 FORMAT(16H GADAP:I TOO BIG)                                       
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   


      SUBROUTINE GADAPs(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAPs1/ NUM,IFU                                            
      EXTERNAL F                                                        
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
    1 FORMAT(16H GADAP:I TOO BIG)                                       
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAP(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAP(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAP(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAP(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAP(X,A0,B0,F,EPS)    
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAP=SUM    
      EPS=EPS/RED   
      RETURN    
      END   


      SUBROUTINE GADAPS2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAPS(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAPS(X,A0,B0,F,EPS)    
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAPS=SUM    
      EPS=EPS/RED   
      RETURN    
      END   

