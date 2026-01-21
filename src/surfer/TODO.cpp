
Scene discussion
===================

- want to use the scene/image rendering with minimal duplicate code in multiple tools:
1. surfer2
2. stitch?
3. a typical solver?
4. lesimage?
5. multiple vars in the same image command

WRITE_IMAGE .... VAR=P,mag(U),T


- debate between more abstraction and registration vs more primative methods with
less awareness of class, but likely more duplicate code to write

- try and remove the vetor math associated with setting up the scene from the WRITE_IMAGE command and the bbox



CtiScene scene;

void processWriteImage(Param * param,const int step) {

  CtiScene scene(param); // WRITE_IMAGE...
  if (scene.write_interval%step != 0)
  return;



  scene.setBbox(xyzmin,xyzmax);

  scene.addTri(x,y,z,...);
  scene.bufferImage();
  scene.writeImage();

  // ......

  scene.setBbox();

  scene.init(param);


  // ALSO, want to reuse parts of imaging for efficient planar rendering, or renering time-dep isosurface in geometry
  // snapshot processing...

  Scene imag1(param);
  imag1.render...


  Scene imag1(param);


  //

  for (step = 0 to 1000) {

    // read new data

    Scene image1_final = image1;
    image1_final.addPlaneData(phi);
    image1_final.writeImage()...





    2. image metadata
    - should all metadat be written to all images? with blanks where nothing is set?
    - how to set metadat from the calling process?
    - range is like the bbox -- may have to come from somewhere?
    - default range is the local range of pixel data


    3. WRITE_IMAGE syntax

    - surface vs volume vars?
    - camera/target vs geom/plane


    e.g.
    WRITE_IMAGE NAME=phi0 GEOM=PLANE x y z nx ny nz VAR=phi,T,mag(u)
    - center the plane on the view and put the plane at depth zero
    - render all geometry behind the plane
    - associate VAR with preceeding GEOM to allow multiple planes, or combinations of planes and surface data?

    WRITE_IMAGE NAME=phi2 BLANK PLANE x y z nx ny nz
    - what view doe sthis one default to?


    WRITE_IMAGE NAME tau_wall GEOM=FAZONE x0,x1,y0 VAR=tau
    - no camera data here, what geometry does this show and what view?

    what about getting rid of GEOM=FAZONE:

    tauwall(surface),pressure(volume) on 3 zones x1,x2,y0 with a default view:

    WRITE_IMAGE NAME=mydata SHOW_ZONE x1,x2,y0 WITH_DATA tauwall,p SHOW_ZONE x3,y3

    WRITE_IMAGE NAME=mydata SHOW_DATA tauwall,p

    - data vs no data
    - surface vs planes vs other volumetric rendering
    - mesh vs no mesh - is mesh a surface var?


    target visualizations:


    PLANE data + non-blanked part of the model
    SURFACE data + other surfaces in gold/gray
    SURFACE data + a plane of data
    SURFACE MESH
    PLANE MESH - stitch
    ISOsurface plus some surfaces in gold/gray

    Tasks/ownership:
    ===========================
    1. ME own SimpleSurface_*.cpp
    - migrate all the junk from surfer2 - e.g. st_flag,grouping,graph of relationships between surfaces,

    2. PHUC - client side architecture, support temporary context menu request
    - client side requests - e.g. command plus write image, not necc, json request again - handshaking

    3. DP - image and scene reduction and working across a few eample tools, normal rendering
    - work with phuc on view control - axis
    - surface data viz
    - data probing?

    4. - help DP with scene/CtiCanvas condensation
    - prototype function registration for processing parameter data
    - interaction with restart file, give access to volume data
    - ME on surface algos














    HIDE_ZONE 13,242,43

    allow rendering of spray data/particles

    flags:
    SHOW_OPEN_EDGES
    SHOW_SURFACE_MESH
    SHOW_SURFACE_VAR
