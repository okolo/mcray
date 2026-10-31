!=============================================================================!
      double precision function spec_int(ep,es,id,reac)
!-----------------------------------------------------------------------------
! interpolation of lab. energy spectra of photons and antiprotons (+antineutrons)
! input:
!   ep  - incident energy per nucleon (GeV);
!   es  - secondary particle energy (GeV);
!   id  - secondary particle type (0 - photon, 1 - antiproton)
!   reac - production process (1 - p+p; 2 - p+He; 3 - He+p)
!   reac=0: photon spectra for p+p collisions based on the non-diffractive part
!       of the Kamae parametrization for ep<ethr and present tables for ep>ethr 
!-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision kamae_nd
      integer reac
      parameter(ngamdat=32,napdat=28,nebin=19,nreac=3)
      dimension xlim(2),ilim(2),wi(3),wj(3)
      common/data/data_gam(ngamdat,nebin,nreac)
     *,data_ap(napdat,nebin,nreac)
      common/ethr/ethr
      data xlim/.55d0,.45d0/
      data ilim/28,23/
      
      if(es.gt.ep.or.es.lt.0.d0)then
       write(*,*)'respect kinematic limits!'
       stop
      endif
      
      if(id*(id-1).ne.0)then
       write(*,*)'no tables for this particle!'
       stop
      endif
      
      if(reac.eq.0.and.ep.lt.ethr)then !for reac=0 and ep<ethr use ND Kamae param.
       if(id.ne.0)stop'reac=0 - for p+p->gamma only!'
       spec_int=kamae_nd(ep,es)
       return
      endif
      
      if(reac*(reac-1)*(reac-2)*(reac-3).ne.0)then
       write(*,*)'no tables for this reaction!'
       stop
      endif
      
      if(ep.lt.10.d0.and.reac>1)then
       write(*,*)'results are questionable below 10 GeV lab.!'
       stop
      endif
      
      if(ep.gt.1.d8)then
       write(*,*)'extrapolation beyond the tabulated range!'
       stop
      endif
      
      if(ep.lt.1.d3)then
       yy=dlog10(ep)*4.d0-3.d0
      else
       yy=dlog10(ep)*2.d0+3.d0
      endif
      jy=max(1,int(yy))
      jy=min(jy,nebin-2)
      if(jy.eq.8)jy=7
      wj(2)=yy-dble(jy)
      wj(3)=wj(2)*(wj(2)-1.d0)*.5d0
      wj(1)=1.d0-wj(2)+wj(3)
      wj(2)=wj(2)-2.d0*wj(3)

      x=es/ep
      if(x.lt.xlim(id+1))then
       xx=50.d0*x+.5d0
      else
       xx=ilim(id+1)+(x-xlim(id+1))*10.d0
      endif
      ix=max(1,int(xx))
      if(ix.eq.ilim(id+1)-1)ix=ix-1
      
      spec_int=0.d0
      ir=max(1,reac)
      if(id.eq.0)then
       ix=min(ix,ngamdat-3)
       if(data_gam(ix,jy,ir)*data_gam(ix,jy+1,ir)
     * *data_gam(ix,jy+2,ir).eq.0.d0)ix=ix-1
      elseif(id.eq.1)then
       ix=min(ix,napdat-3)
       if(data_ap(ix,jy,ir)*data_ap(ix,jy+1,ir)
     * *data_ap(ix,jy+2,ir).eq.0.d0)then
        if(data_ap(ix+1,jy,ir)*data_ap(ix+1,jy+1,ir)
     *  *data_ap(ix+1,jy+2,ir).eq.0.d0)return
        ix=ix+1
       endif
       if(data_ap(ix+2,jy,ir)*data_ap(ix+2,jy+1,ir)
     * *data_ap(ix+2,jy+2,ir).eq.0.d0)ix=ix-1
       if(data_ap(ix+2,jy,ir).eq.0.d0)then
        if(data_ap(ix+1,jy,ir).eq.0.d0)return
        ix=ix-1
       endif
      endif
      wi(2)=xx-dble(ix)
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)
      
      if(id.eq.0)then
       do i=1,3
       do j=1,3
        spec_int=spec_int+log(data_gam(ix+i-1,jy+j-1,ir))*wi(i)*wj(j)
       enddo
       enddo
      elseif(id.eq.1)then
       do i=1,3
       do j=1,3
        if(ix.eq.ilim(id+1)-1.and.i.eq.1)then
         dspec=log(data_ap(ix+i-5,jy+j-1,ir))
        else
         dspec=log(data_ap(ix+i-1,jy+j-1,ir))
        endif
        spec_int=spec_int+dspec*wi(i)*wj(j)
       enddo
       enddo
      endif
      spec_int=exp(spec_int)
      return
      end
      
!=============================================================================!
      subroutine spec_ini
!-----------------------------------------------------------------------------
! read tables
!-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(ngamdat=32,napdat=28,nebin=19,nreac=3)
      common/data/data_gam(ngamdat,nebin,nreac)
     *,data_ap(napdat,nebin,nreac)
      common/ethr/ethr

      ethr=5.d1 !threshold for switching to the ND Kamae param. (for reac=0 only)
  
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!                                                   !'
      write(*,*) '! ppfrag -- for more infor see                      !'
      write(*,*) '! http://sourceforge.net/projects/ppfrag            !'
      write(*,*) '!                                                   !'
      write(*,*) '! if you use data from reac=1,2,3, please cite      !'
      write(*,*) '! M.Kachelriess, S.Ostapchenko                      !'
      write(*,*) '! arxiv[1206.xxxx], submitted to PRD                !'
      write(*,*) '!                                                   !'
      write(*,*) '! if you use reac=0, please cite additionally       !'
      write(*,*) '! T. Kamae et al., ApJ 647 (2006) 692.              !'
      write(*,*) '!                                                   !'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         
      open(1,file='gamfrag.dat',status='old')
      read(1,*)data_gam
      close (1)
      open(1,file='apfrag.dat',status='old')
      read(1,*)data_ap
      close (1)
      write(*,*)'cross section tables read'
      return
      end
      
!=============================================================================!
      double precision function kamae_nd(ep,es)
!-----------------------------------------------------------------------------
! parametrization of lab. energy spectra of photons in non-diffractive 
! proton-proton interactions from T. Kamae et al., ApJ 647 (2006) 692.
! input:
!   ep  - incident energy per nucleon (GeV);
!   es  - secondary photon energy (GeV)
!-----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      
      amn=.938d0
      ekin=ep-amn
      y=dlog10(ekin)-3.d0      
      if(ekin.gt..487d0)then
       a0g=-0.51187d0*(y+3.3d0)+7.6179d0*(y+3.3d0)**2
     * -2.1332d0*(y+3.3d0)**3+0.22184d0*(y+3.3d0)**4
       a1g=-1.2592d-5+1.4439d-5*exp(-0.29360d0*(y+3.4d0))
     * +5.9363d-5/(y+4.1485d0)+2.2640d-6*y-3.3723d-7*y**2
       a2g=-174.83d0+152.78d0*dlog10(1.5682d0*(y+3.4d0))
     * -808.74d0/(y+4.6157d0)
       a3g=0.81177d0+0.56385d0*y+0.0040031d0*y**2-0.0057658d0*y**3
     * +0.00012057d0*y**4
       a4g=0.68631d0*(y+3.32d0)+10.145d0*(y+3.32d0)**2-4.6176d0
     * *(y+3.32d0)**3+0.86824d0*(y+3.32d0)**4-0.053741d0*(y+3.32d0)**5
       a5g=9.0466d-7+1.4539d-6*dlog10(0.015204d0*(y+3.4d0))
     * +1.3253d-4/(y+4.7171d0)**2-4.1228d-7*y+2.2036d-7*y**2
       a6g=-339.45d0+618.73d0*dlog10(0.31595d0*(y+3.9d0))
     * +250.20d0/(y+4.4395d0)**2
       a7g=-35.105d0+36.167d0*y-9.3575d0*y**2+0.33717d0*y**3
       a8g=0.17554d0+0.37300d0*y-0.014938d0*y**2+0.0032314d0*y**3
     * +0.0025579d0*y**4
      endif
     
      if(ekin.lt.1.95d0)then
       ry=3.05d0*exp(-107.d0*((y+3.25d0)/(1.d0+8.08d0*(y+3.25d0)))**2)
      else
       ry=1.01d0
      endif

      if(ekin.gt..488d0.and.ekin.lt.1.95d0)then
       c0g=2.4316d0*exp(-69.484d0*((y+3.1301d0)/(1.0d0+0.14921d0
     * *(y+3.1301d0)))**2)-6.3003-9.5349d0/y+0.38121d0*y**2
       c1g=56.872d0+40.627d0*y+7.7528d0*y**2
       c2g=-5.4918d0-6.7872d0*tanh(4.7128d0*(y+2.1d0))+0.68048d0*y
       c3g=-0.36414d0+0.039777d0*y
       c4g=-0.72807d0-0.48828d0*y-0.092876d0*y**2
      endif
      
      if(ekin.gt..69d0.and.ekin.lt.2.76d0)then
       d0g=3.2433d0*exp(-57.133d0*((y+2.9507d0)/(1.0d0+1.2912d0
     * *(y+2.9507d0)))**2)-1.0640d0-0.43925d0*y
       d1g=16.901d0+5.9539d0*y-2.1257d0*y**2-0.92057d0*y**3
       d2g=-6.6638d0-7.5010d0*tanh(30.322d0*(y+2.1d0))+0.54662d0*y
       d3g=-1.50648d0-0.87211d0*y-0.17097d0*y**2
       d4g=0.42795d0+0.55136d0*y+0.20707d0*y**2+0.027552d0*y**3
      endif

      x=dlog10(es)
      if(ekin.gt..487d0)then
       fnd=a0g*exp(-a1g*(x-a3g+a2g*(x-a3g)**2)**2)
     * +a4g*exp(-a5g*(x-a8g+a6g*(x-a8g)**2+a7g*(x-a8g)**3)**2)
      else
       fnd=0.d0
      endif

      if(ekin.gt..488d0.and.ekin.lt.1.95d0)then
       fdelt=c0g*exp(-c1g*((x-c2g)
     * /(1.d0+c3g*(x-c2g)+c4g*(x-c2g)**2))**2)
      else
       fdelt=0.d0
      endif
      
      if(ekin.gt..69d0.and.ekin.lt.2.76d0)then
       fres=d0g*exp(-d1g*((x-d2g)
     * /(1.d0+d3g*(x-d2g)+d4g*(x-d2g)**2))**2)
      else
       fres=0.d0
      endif
      
      kamae_nd=max(0.d0,fnd*ry/(1.d0+exp(15.d0*(-2.6d0-x)))
     */(1.d0+exp(44.d0*(x-.96d0*dlog10(ekin)))))
      kamae_nd=kamae_nd+max(0.d0,fdelt
     */(1.d0+exp(75.d0*(x-dlog10(ekin)))))
      kamae_nd=kamae_nd+max(0.d0,fres
     */(1.d0+exp(75.d0*(x-dlog10(ekin)))))
      return
      end
