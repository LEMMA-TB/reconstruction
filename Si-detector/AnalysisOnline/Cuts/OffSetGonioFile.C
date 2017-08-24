///////////////////////////////////////////
//          OFFSET GONIO FILES           //
//      for Coherent data analysis       //
//     01 July 2012 - Enrico Bagli       //
///////////////////////////////////////////



/////////******* ST2012_1 ********///////// - EB
if(runnumb==97412) offsetX = (9.2);
/////////*************************/////////



/////////****** GE_UNICOR_2 ******///////// - EB
if((runnumb==97422) || (runnumb==97423)) offsetX = (+85.);
if(runnumb==97426) offsetX = (32.);
if(runnumb==97427) offsetX = (2.034);
/////////*************************/////////



/////////********** SiGe1 ********///////// - EB
if((runnumb==97497) || (runnumb==97500)) offsetX = (-2.443 + 0.6726 + gonio_rot-3180330. +2.88 +1.78 -0.06);
/////////*************************/////////



/////////********* ST27 **********///////// - EB
if((runnumb>=97660 && runnumb<=97662 && N_RUNS>2) || (runnumb==97684 && N_RUNS>2))
{
    offsetX = gonio_rot - 3183783.507;
    offsetY = gonio_crad - 292595.649;
}
if(runnumb>=97666 && runnumb<=97681 && N_RUNS>2)
{
    offsetX = gonio_rot - 3183778.18;
    offsetY = gonio_crad - 294171.;
}
if( (runnumb==97684 && N_RUNS==1) || ( (runnumb==97684 || runnumb == 97691) && N_RUNS==2 ) )
{
    offsetX = - 0.922 + gonio_rot - 0.31837832E+07;
}
/////////*************************/////////



/////////******* ST2012_1 ********///////// - EB
if(runnumb>=97703 && runnumb<=97706)
{
    offsetX = -(gonio_rot - 3176955.);
    offsetY = gonio_crad - 292672.;
}

if(runnumb>=97713 && runnumb<=97715)
{
    offsetX = - (gonio_rot - 3179780.);
}

if(runnumb>=97721 && runnumb<=97737 && N_RUNS>2)
{
    offsetX = - (gonio_rot - 3176959.);
    offsetY = gonio_crad - 291017.;
}
/////////*************************/////////



/////////******* ST2012_2 ********///////// - EB
if(runnumb>=97742 && runnumb<=97743 && CENTER && !GONIO_ZERO)
{
    offsetX = gonio_rot + 39.94;
    offsetY = gonio_crad + 39.94;
}
if(runnumb == 97744) offsetX = -6.9;
/////////*************************/////////



/////////********** SiGe2 ********///////// - EB
if(runnumb==97785) offsetX = (+1.59);
/////////*************************/////////



/////////****** GE_UNICOR_2 ******///////// - EB
if(runnumb == 97802) offsetX = +101.5;
if(runnumb==97803)  offsetX =  -3.75	;
/////////*************************/////////

/////////******* STGIGI_01 *******///////// - EB
// PIANO 111
// ANGOLO 200 urad
// 2 x 2 x 55
// TORSION 0
if(runnumb==100443 ||
   runnumb==100444 ||
   runnumb==100445 ||
   runnumb==100446 ||
   runnumb==100448){
    offsetX = -21. - (gonio_rot + 308468.);
}
//FIXED GONIO CH POSITION
if(runnumb==100450
   || runnumb==100455
   || runnumb==100457
   || runnumb==100459
   || runnumb==100461
   || runnumb==100463
   ){
    //offsetX = (gonio_rot - 308468.1);
  offsetX =  (-21.28)  +5.176763 ;
  //thXin *= -1.;
  //thXout *= -1.;
}
if(runnumb==100475){
    offsetX =  -21.;
}
//MOVED GONIO +80.
if(runnumb>=100476 &&
   runnumb<=100506){
    //offsetX =  +59.+6.5-3.14896+5.9338;
}
//MOVED GONIO +75.
if(runnumb==100507){
    offsetX =  +59.+6.5+75.;
}
//MOVED GONIO -150.
if(runnumb>=100508 &&
   runnumb<=100514){
    offsetX =  -21.28+(308468.1-308463.) +5.176763;
}
if(runnumb==100520){
    offsetX =  -50.-3.62;
}
//MOVED GONIO +30
if(runnumb>=100521 &&
   runnumb<=100532){
    offsetX =  -20.-3.62;
}

//INSTALLED QM
if(runnumb>=100534){
    offsetX =  0.;
}

if(runnumb==100537){
    offsetX =  90.;
}

if(runnumb == 90127) {
    offsetX = -49.2422;
}


/////////*************************/////////

if(runnumb == 1724) {
    offsetX = -9.76;
}

if(runnumb == 2726){
    offsetX = -1.5088 + 1.-1.;
}

if(runnumb == 3078){
    offsetX = -9.097/2.-0.9;
}

if(runnumb == 2753){
    offsetX = +1.35;
}

/////////*************************/////////

if(runnumb == 110138 || runnumb == 110138 ){
    offsetX = +20.-15.75-2.28;
}

if(runnumb == 110161 ){
    offsetX = +45.-6.5 - 27.68 + 0.5;
}

if(runnumb == 110170 ||
   runnumb == 110172 ||
   runnumb == 110174 ||
   runnumb == 110176 ||
   runnumb == 110178){
    offsetX = +45.-6.5 - 27.68 + 0.5;
//    offsetX = +45.-6.5-1.6;
}

if(runnumb == 112359 ||
   runnumb == 112363){
    offsetX = -40. -28.;
}

if(runnumb >= 112392 && runnumb <= 112413){
    offsetX = 7. + 5.;
    offsetY = 133.;
}

if(runnumb >= 112432 && runnumb <= 112437){
    offsetX = 7. + 5.;
    offsetY = 133.;
}

if(runnumb >= 112438 && runnumb <=112454) {
    offsetX = 7. + 5. + 50.;
    offsetY = 133. - 115.;
}

if(runnumb >= 112455){
    offsetX = 7. + 5.;
    offsetY = 133.;
}

if(runnumb >= 112460 && runnumb <=112490) {
    offsetX = 0.;// 7. + 5. + 50.;
    offsetY = 0.;// 133. - 115.;
}
if(runnumb >= 112491 && runnumb <=112520) {
    offsetX = 0.;// 7. + 5. + 50.;
    offsetY = 0.;// 133. - 115.;
}


