{ MAGNET_COIL.PDE  
  

 
  According to Maxwell's equations, 
        curl H = J 
        div B = 0 
        B = mu*H or H = B/mu
 
  where B is the manetic flux density
        H is the magnetic field strength
        J is the electric current density
  and  mu is the magnetic permeability of the material.
 
  The magnetic vector potential A is related to B by 
        B = curl A 
  therefore 
        curl( (1/mu)*curl A ) = J
    if rmu = 1/mu
    then 
          curl( (rmu)*curl A ) = J    
 
  This equation is usually supplmented with the Coulomb Gauge condition 
        div A = 0.
 
  In the axisymmetric case, the current is assumed to flow only in the
  azimuthal direction, and only the azimuthal component of the vector
  potential is present.  Henceforth, we will simply refer to this component as A.
 
  The Coulomb Gauge is identically satisfied, and the PDE to be solved in this
  model takes the form 
        curl((1/mu)*curl (A)) = J(x,y)     in the domain
                           A  = g(x,y)     on the boundary.
 
  The magnetic induction B takes the simple form 
        B = (-dz(A), 0, dr(A)+A/r)
 
  and the magnetic field is given by 
        H = (-dz(A)/mu, 0, (dr(A)+A/r)/mu)
 
  Expanding the equation in cylindrical geometry results in the final equation, 
        dz(dz(A)/mu) + dr((dr(A)+A/r)/mu) = -J
 
  The interpretation of the natural boundary condition becomes 
        Natural(A) = n X H
 
  where n is the outward surface-normal unit vector.
 
  Across boundaries between regions of different material properties, the
  continuity of (n X H) assumed by the Galerkin solver implies that the
  tangential component of H is continuous, as required by the physics.
 
 
  In this simple test problem, we consider a circular coil whose axis of
  rotation lies along the X-axis. We bound the coil by a distant spherical
  surface at which we specify a boundary condition (n X H) = 0.
  At the axis, we use a Dirichlet boundary condition A=0.
 
}
 
title 'AXI-SYMMETRIC MAGNETIC FIELD'
 
coordinates
    { Cylindrical coordinates, with cylinder axis along Cartesian X direction }
    xcylinder(Z,R)  
 
variables
    Aphi        { the azimuthal component of the  vector potential }
 
definitions
   Lx =0.1
   Ly = 0.1
    mu0 = pi *4e-7            { the permeability }
    mu 
    J              { the source defaults to zero }
    rmu = 1/ mu
    nudge = 0.0001
    wire_length = 0.02
    pen_len = 2/100
    
    a1 = 0.0179		                                                                                        ! Inner radius
    a2 = 0.0214                                                                                                  ! outer Radius
    sth =0.0015                                                                                                   ! Thickness of the skin
    ar = 0.05                                                                                                        ! radius of the arm
    fth = 0.0085                                                                                                   ! Thickness of the Fat 
    mth=0.0275                                                                                                   ! Thickness of muscle 
    ! Radius  of the wire
    muu_skin = pi *4e-7 
    muu_muscle=pi *4e-7 
    muu_fat = pi *4e-7 
    muu_wire = pi *4e-7 
    
     pvc_s_d = (8.44/1000)/2
     wire_rad = pvc_s_d
     cop_s_d = (6.3/1000)/2
     sil_s_d = (5/1000)/2
     di_s_d = (4.7/1000)/2
     cop_c_d = (0.74/1000)/2
  
    current = 2.14
    current_density = current/(pi*wire_rad^2)
    Bz = dr(r*Aphi)/r                                                                                       ! The Z-component of the Flux density
    Br = dz(Aphi)                                                                                            ! The R-component of the Flux density Note: the phi-component
    Bmag = sqrt(Bz^2 + Br^2)                                                                         ! The Magnitude of the flux density
    H = curl(Aphi)/mu                                                                                      ! The magnetic field strength
    Hmag = magnitude(H)                                                                               ! The magnitude magnetic field strength
    B = curl(Aphi)                                                                                            ! The magnetic flux density
    !Bmag = magnitude(B)
    !J1 = curl(H)
initial values
    Aphi = 2            { unimportant unless mu varies with H }
 
equations
    { FlexPDE expands CURL in proper coordinates }
    Aphi : curl(rmu*curl(Aphi)) = J 
 
boundaries
    region 1
    J = 0
     mu = mu0
       START 'ring' ( -Lx , 0)
           line to (0,Lx)
           line to (Ly, Lx)
           line to (0, Ly)
           line to close 
 
    region 2
        mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  wire_rad)
           line to close 
    {
    region 3
        J =0 
    	mu = muu_skin
        start 'Skin'((a2-a1)/2,ar +nudge)
        line to ((a2-a1)/2,  wire_rad ) to ((a2-a1)/2+sth,  wire_rad ) to ((a2-a1)/2+sth,ar+ nudge) to close
        }
   region 4
        J  = 0
    	mu = muu_fat
        start 'Fat'( (a2-a1)/2  +sth ,ar + nudge)
        line to ((a2-a1)/2 + sth , wire_rad) to ((a2-a1)/2+ sth +fth, wire_rad) to ((a2-a1)/2 + sth +fth,ar+ nudge) to close
    region 5
         J = 0
    	mu = muu_muscle 
        start 'Muscle'( (a2-a1)/2  +sth + fth ,ar + nudge)
        line to ((a2-a1)/2 + sth + fth , wire_rad ) 
        line to ((a2-a1)/2+ sth +fth+ mth,  wire_rad ) 
        line to ((a2-a1)/2 + sth +fth +mth,ar + nudge) 
        line to close
        
          start 'Muscle'( pen_len+wire_rad, nudge)
        line to (pen_len+wire_rad,  wire_rad ) 
        line to ((a2-a1)/2+ sth +fth+ mth,   wire_rad ) 
        line to ((a2-a1)/2 + sth +fth +mth, nudge) 
        line to close
        
    region 6
        mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, cop_s_d )
           line to (-wire_length , cop_s_d)
           line to (-wire_length ,  wire_rad)
           line to close 
   { 
    region 7
               mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  wire_rad)
           line to close 
    region 8
            mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  wire_rad)
           line to close 
    
    region 9
            mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  wire_rad)
           line to close 
    
    region 10
            mu = muu_wire
        J = current_density
       START 'ring' ( wire_rad +pen_len , wire_rad )
           line to ( wire_rad + pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  wire_rad)
           line to close 
    }
monitors
    contour(Bz) zoom(-2,0,4,4) as 'FLUX DENSITY B'
    contour(Aphi) as 'Potential'
 
plots
    grid(z,r)
    contour(Bz)  as 'FLUX DENSITY B'
    contour(Bz) zoom(-2,0,4,4)  as 'FLUX DENSITY B'
    contour(Bmag) as 'Magnitude of the Flux Density (T)' PNG
    contour(Hmag) as  'Magnitude of Field intensity (A/m)' PNG
    elevation(Aphi,dr(Aphi),Aphi/r,dr(Aphi)+Aphi/r,Aphi+r*dr(Aphi)) 
        from (0,0) to (0,1) as 'Bz'
    vector(dr(Aphi)+Aphi/r,-dz(Aphi)) as 'FLUX DENSITY B'
    vector(dr(Aphi)+Aphi/r,-dz(Aphi)) zoom(-2,0,4,4) as 'FLUX DENSITY B'
    vector(H) as 'Vector Plot of magnetic field intensity H (A/m)' PNG
    vector(B) as 'Vector plot of the Manetic Flux Density B (T)' PNG
    vector(Aphi)  as 'VECTOR PLOT OF MAGNETIC POTENTIAL' PNG
    contour(Aphi)  as 'MAGNETIC POTENTIAL' PNG
    contour(Aphi) zoom(-2,0,4,4)  as 'MAGNETIC POTENTIAL'
    
    surface(Aphi)  as 'MAGNETIC POTENTIAL'  viewpoint (-1,1,30)
    
Summary
    !report (GlobalMax(Bmag, 1)) as 'Maximum flux density  (Wb/m^2)'
    !report (GlobalMax(Bmag, 2)) as 'Maximum flux density in the wire (Wb/m^2)'
    report (GlobalMax(Bmag, 3)) as 'Maximum flux density in the skin (Wb/m^2)'
    report (GlobalMax(Bmag, 4)) as 'Maximum flux density in the fat (Wb/m^2)'
    report (GlobalMax(Bmag, 5)) as 'Maximum flux density in the muscle (Wb/m^2)'
    
    !report (GlobalMax(Hmag, 1)) as 'Maximum Magnetic Field intensity (A/m)'
    !report (GlobalMax(Hmag, 2)) as 'Maximum Magnetic Field intensity in the wire (A/m)'
    report (GlobalMax(Hmag, 3)) as 'Maximum Magnetic Field intensity in the skin (A/m)'
    report (GlobalMax(Hmag, 4)) as 'Maximum Magnetic Field intensity in the fat (A/m)'
    report (GlobalMax(Hmag, 5)) as 'Maximum Magnetic Field intensity in the muscle (A/m)'
end "CIEPcbP1p8UrD0Ykk9rG4yeq+J6b384kTZRnemmylDzh4NSeLTdnTMCRkOBMoS23tpclJAzvcbEpb3QAb6K7Qovjmwqs32XquSq7Zn/zuZvD0C8taSgsHPSQzqwzfwJ8wpGqtCfyHAV+cyJ02PPcvZ8Sc5HVbAp95qxUD3S3CqH"
