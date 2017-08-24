///////////////////////////////////////////
//               CUT FILES               //
//      for Coherent data analysis       //
//     01 July 2012 - Enrico Bagli       //
///////////////////////////////////////////

double torsionGiGi1 = 0.99-2.786;


/////////******* ST2012_1 ********///////// - EB
if((runnumb >= 97405) &&
   (runnumb <= 97412)){
    deflXmar1=-64.;    deflXmar2=+512.;
    deflYmar1=-64.;    deflYmar2=+512.;
    torsion_posYthX=10.;//murad on mm
    X_CUT_left = 10500.;    X_CUT_right = 11300.;
    Y_CUT_left = 8500.+1500.;    Y_CUT_right = 12500.+1500.;
}
/////////*************************/////////



/////////****** GE_UNICOR_2 ******///////// - EB
else if((runnumb >= 97413) &&
        (runnumb <= 97429)){
    deflXmar1=-128.;    deflXmar2=+384.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=-2.8;//murad on mm
    X_CUT_left = 10610.-300.;    X_CUT_right = 10610.+300.;
    Y_CUT_left = 11500.;    Y_CUT_right = 12500.;
}
/////////*************************/////////



/////////********** MSTN1 ********///////// - EB
else if((runnumb >= 97430) &&
        (runnumb <= 97475)){
    deflXmar1=-256.;    deflXmar2=+384.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=16.;//murad on mm
    X_CUT_left = 11050.;    X_CUT_right = 11300.;
    Y_CUT_left = 11500.;    Y_CUT_right = 12500.;
}
else if((runnumb >= 97476) &&
        (runnumb <= 97482)){
    deflXmar1=-256.;    deflXmar2=+384.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=-16.;//murad on mm
    X_CUT_left = 10400.;    X_CUT_right = 10500.;
    Y_CUT_left = 11500.;    Y_CUT_right = 12500.;
}
/////////*************************/////////



/////////********** SiGe1 ********///////// - EB
else if((runnumb >= 97483) &&
        (runnumb <= 97500)){
    
    deflXmar1=-256.;    deflXmar2=+384.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=-16.;//murad on mm
    X_CUT_left = 10807.-235.;    X_CUT_right = 10807.+235.;
    Y_CUT_left = 11500.-250.;    Y_CUT_right = 11500.+250.;
}
/////////*************************/////////



/////////********** MSTN2 ********///////// - EB
else if((runnumb >= 97594) &&
        (runnumb <= 97601)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-16.;//murad on mm
    X_CUT_left = 10650.;    X_CUT_right = 10975.;
    Y_CUT_left = 11500.;    Y_CUT_right = 12500.;
}
/////////*************************/////////



/////////******* MST-SPS *********///////// - EB
else if((runnumb >= 97602) &&
        (runnumb <= 97621)){
    deflXmar1=-1024.;    deflXmar2=+256.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-16.;//murad on mm
    X_CUT_left = 10635.;    X_CUT_right = 10835.;
    Y_CUT_left = 11500;    Y_CUT_right = 12500;
}
/////////*************************/////////



/////////**** GE_UNICOR_DISLO ****///////// - EB
else if((runnumb >= 97622) &&
        (runnumb <= 97634)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-16.;//murad on mm
    X_CUT_left = 10430.;    X_CUT_right = 10620.;
    Y_CUT_left = 11500;    Y_CUT_right = 12500;
}
/////////*************************/////////



/////////******* UNDULATOR *******///////// - LB
else if((runnumb >= 97635) &&
        (runnumb <= 97653)){
    deflXmar1=-512.;    deflXmar2=+512.;
    deflYmar1=-512.;    deflYmar2=+512.;
    torsion_posYthX=0.;//murad on mm
    X_CUT_left = 10515;    X_CUT_right = 10550;
    Y_CUT_left = 11780-1000;    Y_CUT_right = 11780 + 1000;
    //Y_CUT_left = 10500;    Y_CUT_right = 13000;
}
/////////*************************/////////



/////////********* ST27 **********///////// - EB
else if((runnumb >= 97654) &&
        (runnumb <= 97691)){
    deflXmar1=-128.;    deflXmar2=+128.;
    deflYmar1=-128.;    deflYmar2=+128.;
    torsion_posYthX = -9.99;//murad on mm
    torsion_posXthY = -160.;//murad on mm
    X_CUT_left = 10595.-400.;    X_CUT_right = 10595.+400.;
    Y_CUT_left = 11500;    Y_CUT_right = 12500;
}
/////////*************************/////////



/////////******* ST2012_1 ********///////// - EB
else if((runnumb >= 97692) &&
        (runnumb <= 97702)){
    
    deflXmar1=-128.;    deflXmar2=+512.;
    deflYmar1=-128.;    deflYmar2=+128.;
    torsion_posYthX=20;//murad on mm
    torsion_posXthY = +260.;
    X_CUT_left = 10200.;    X_CUT_right = 11000.;
    Y_CUT_left = 11000;    Y_CUT_right = 13000;
}
else if((runnumb >= 97703) &&
        (runnumb <= 97720))
{
    thXin = 64;
    thYin = 64;
    deflXmar1=-256.;    deflXmar2=+512.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=0.5;//murad on mm
    torsion_posXthY=+260.;
    X_CUT_left = 10370.-400.;    X_CUT_right = 10370.+400.;
    Y_CUT_left = 11000;    Y_CUT_right = 13000.;
}
else if((runnumb >= 97721) &&
        (runnumb <= 97737))
{
    thXin = 128;
    thYin = 128;
    deflXmar1=-256.;    deflXmar2=+512.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=0.5;//murad on mm
    torsion_posXthY=+260.;
    X_CUT_left = 10370.-400.;    X_CUT_right = 10370.+400.;
    Y_CUT_left = 11000;    Y_CUT_right = 13000.;
}
/////////*************************/////////



/////////******* ST2012_2 ********///////// - EB
else if((runnumb >= 97738) &&
        (runnumb <= 97746)){
    
    deflXmar1=-128.;    deflXmar2=+384.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-7.446;//murad on mm
    X_CUT_left = 10640.-400.;    X_CUT_right = 10640.+400.;
    Y_CUT_left = 12000.-1000.;    Y_CUT_right = 12000.+1000.;
}
/////////*************************/////////



/////////******* ST2012_3 ********///////// - EB
else if((runnumb >= 97749) &&
        (runnumb <= 97752)){
    
    deflXmar1=-128.;    deflXmar2=+896.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX = -9;//murad on mm
    X_CUT_left = 10230;    X_CUT_right = 11050;
    Y_CUT_left = 11000;    Y_CUT_right = 13000;
}
else if((runnumb >= 97753) &&
        (runnumb <= 97768)){
    
    deflXmar1=-128.;    deflXmar2=+896.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX = 2.7;//murad on mm
    X_CUT_left = 10230;    X_CUT_right = 11050;
    Y_CUT_left = 11000;    Y_CUT_right = 13000;
}
/////////*************************/////////



/////////********** SiGe2 ********///////// - EB
else if((runnumb >= 97769) &&
        (runnumb <= 97787)){
    
    deflXmar1=-400.;    deflXmar2=+400.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=2.7;//murad on mm
    X_CUT_left = 10723. -235. ;    X_CUT_right = 10723. +235. ;
    Y_CUT_left = 11500. - 250.;    Y_CUT_right = 11500. + 250.;
    
}
/////////*************************/////////



/////////****** GE_UNICOR_2 ******///////// - EB
else if((runnumb >= 97789) &&
        (runnumb <= 97807)){
    
    deflXmar1=-256.;    deflXmar2=+1024.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-3.;//murad on mm
    X_CUT_left = 10450;    X_CUT_right = 10900;
    Y_CUT_left = 11000;    Y_CUT_right = 13000;
}
/////////*************************/////////



/////////******* ST2012_4 ********///////// - EB
else if((runnumb >= 97808) &&
        (runnumb <= 97999)){
    
    deflXmar1=-256.;    deflXmar2=+1024.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=-12.5;//murad on mm
    X_CUT_left = 10250-288;    X_CUT_right = 11100-288;
    Y_CUT_left = 11000;    Y_CUT_right = 13000;
}
/////////*************************/////////


/////////* DECHANNELING ANALYSIS *///////// - EB
//if(runnumb==97684) divAcutX = 7.035;
//if(runnumb==97744) divAcutX = 6.671;
//if(runnumb==97717) divAcutX = 6.464;
//if(runnumb==97767) divAcutX = 5.916;
//if(runnumb==97825) divAcutX = 5.027;
/////////*************************/////////




/////////******* STGIGI_01 *******///////// - EB
// PIANO 111
// ANGOLO 200 urad
// 2 x 2 x 55
// TORSION 0

// ELECTRON FROM 100402 TO 100463


else if((runnumb >= 100403) &&
        (runnumb <= 100435)){
    deflXmar1=-256.;    deflXmar2=+1024.;
    deflYmar1=-256.;    deflYmar2=+256.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8412.-50.;    X_CUT_right = 9262.+50.;
    Y_CUT_left = 0;    Y_CUT_right = 20000;
    
}

else if((runnumb >= 100436) &&
        (runnumb <= 100463)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    //torsion_posYthX=-5.+3.48;//murad on mm
    //X_CUT_left = 8412.-50.+180.;    X_CUT_right = 9262.+50.+180.;//during online data analysis
    X_CUT_left = 8970.-944.;    X_CUT_right = 8970.+944.;
    
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=64.;// murad
    thYin=64.;// murad
    meanchang = +140.;
}

// POSITRON FROM 100464 TO 100501

/////////*************************/////////

else if((runnumb >= 100464) &&
        (runnumb <= 100478)){
    deflXmar1=-256.;    deflXmar2=+128.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8412.-50.+180.;    X_CUT_right = 9262.+50.+180.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +140.;
}

//MOVED LINEAR STAGE +250 um

else if((runnumb >= 100479) &&
        (runnumb <= 100501)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    torsion_posXthX=0.;//murad on mm
    X_CUT_left = 9202.-944.;    X_CUT_right = 9202+944.;
    Y_CUT_left = 12790-4000;
    Y_CUT_right = 12790+4000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +140.;
}

// ELECTRON FROM 100502 TO 100

else if((runnumb >= 100502) &&
        (runnumb <= 100505)){
    deflXmar1=-256.;    deflXmar2=+128.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8412.-50.+180.+250.;    X_CUT_right = 9262.+50.+180.+250.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +160.;
}

//MOVED LINEAR STAGE -900 um

else if((runnumb >= 100506) &&
        (runnumb <= 100506)){
    deflXmar1=-256.;    deflXmar2=+128.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8412.-50.+180.+250.-900.;    X_CUT_right = 9262.+50.+180.+250.-900.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +160.;
}

//MOVED LINEAR STAGE +0.07 um
//MOVED GONIO STAGE +75 urad

else if((runnumb >= 100507) &&
        (runnumb <= 100507)){
    deflXmar1=-256.;    deflXmar2=+128.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8412.-50.+180.+250.-900.+75.;    X_CUT_right = 9262.+50.+180.+250.-900.+75.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +160.;
}

//MOVED GONIO STAGE -150 urad

else if((runnumb >= 100508) &&
        (runnumb <= 100513)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=torsionGiGi1;//murad on mm
    X_CUT_left = 8377.-944.;    X_CUT_right = 8377.+944.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = +160.;
}

//////////////////////////////////////////
/////////////////////////////////////////
////////////////////////////////////////
///////////////////////////////////////
//INSTALLED NEW CRYSTAL QM MAZZOLARI SI 01

else if((runnumb >= 100515) &&
        (runnumb <= 100532)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    torsion_posYthX=5.-4.;//murad on mm
    torsion_posXthX=-40.5-2.;//murad on mm
    X_CUT_left = 8412.-50.+180.+250.-900.+75.-2000.;    X_CUT_right = 9262.+50.+180.+250.-900.+75.+2000.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = -130.;
}

//INSTALLED NEW CRYSTAL QM DDS GE 01

else if((runnumb >= 100534) &&
        (runnumb <= 100537)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    //torsion_posYthX=-60.+13.5;//murad on mm
    torsion_posYthX=-62.;//murad on mm
    torsion_posXthX=-14.;//murad on mm
    X_CUT_left = 5500.;    X_CUT_right = 12500.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = -140.;
}

//ELECTRONS

else if((runnumb >= 100538) &&
        (runnumb <= 100900)){
    deflXmar1=-256.;    deflXmar2=+256.;
    deflYmar1=-64.;    deflYmar2=+64.;
    //torsion_posYthX=-60.+13.5;//murad on mm
    torsion_posYthX=-62.;//murad on mm
    torsion_posXthX=-14.;//murad on mm
    X_CUT_left = 8412.-50.+180.+250.-900.+75.-2000.;    X_CUT_right = 9262.+50.+180.+250.-900.+75.+2000.;
    Y_CUT_left = 12790-20000;
    Y_CUT_right = 12790+20000;
    thXin=256.;// murad
    thYin=128.;// murad
    meanchang = -140.;
}

/////////////////////////////////////
/////////  STRIP LiNbO3 y-plane /////
/////////////////////////////////////
if( (runnumb >= 90112) &&
   (runnumb <= 90119) )
{
    deflXmar1=-100.;    deflXmar2=412.;
    deflYmar1=-128.;    deflYmar2=128.;

    thXin=64.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-77.9;//murad on mm
    meanchang = +173.;

    Y_CUT_left = 7600.-1024.;    Y_CUT_right = 7600.+1024.;
}
if( (runnumb >= 90120) &&
   (runnumb <= 90148) &&
   (runnumb != 90126) &&
   (runnumb != 90127))
{
    thXin=256.;    thXout=256.;
    deflXmar1=-100.;    deflXmar2=412.;
    deflYmar1=-128.;    deflYmar2=128.;
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-177.9;//murad on mm
    meanchang = +190.;

    X_CUT_left = 13887.21-411.;    X_CUT_right = 13887.21+411.;
    Y_CUT_left = 7600.-1024.;    Y_CUT_right = 7600.+1024.;
}
if( (runnumb == 90126) || (runnumb==90127))
{
    thXin=256.;    thXout=256.;
    deflXmar1=-100.;    deflXmar2=412.;
    deflYmar1=-128.;    deflYmar2=128.;
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-177.9;//murad on mm
    meanchang = +170.;
    
    X_CUT_left = 13887.21-130.;    X_CUT_right = 13887.21+70.;
    X_CUT_left = 13887.21-411.;    X_CUT_right = 13887.21+411.;
    //X_CUT_left = 13887.21;    X_CUT_right = 13887.21+10.;
    Y_CUT_left = 7350.-512.;    Y_CUT_right = 7350.+512.;
}
if( ((runnumb >= 90131) &&
     (runnumb <= 90132)) ||
   ((runnumb >= 90136) &&
    (runnumb <= 90137)) ||
   ((runnumb == 90139)) ||
   ((runnumb == 90147))  ||
   ((runnumb == 90148)) )
{
    thXin=128.;    thXout=256.;
    deflXmar1=-100.;    deflXmar2=412.;
    deflYmar1=-128.;    deflYmar2=128.;
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-177.9;//murad on mm
    meanchang = +190.;

    X_CUT_left = 13480.;    X_CUT_right = 14220.;
    Y_CUT_left = 7600.-1024.;    Y_CUT_right = 7600.+1024.;
}

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

/////////** UA9 2014 - December **/////////

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

if( (runnumb >= 1634) &&
   (runnumb <= 1637) )
{
    deflXmar1=-800.;    deflXmar2=400.;
    deflYmar1=-400.;    deflYmar2=400.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-77.9;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = -670.;    X_CUT_right = 200.;
    Y_CUT_left = -2048.;    Y_CUT_right = 2048.;
}

if( (runnumb >= 1667) &&
   (runnumb <= 1668) )
{
    deflXmar1=-800.;    deflXmar2=400.;
    deflYmar1=-400.;    deflYmar2=400.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-77.9;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = -750.;    X_CUT_right = -690.;
    Y_CUT_left = -2048.;    Y_CUT_right = 2048.;
}

if( (runnumb >= 1675) &&
   (runnumb <= 1676) )
{
    deflXmar1=-800.;    deflXmar2=400.;
    deflYmar1=-400.;    deflYmar2=400.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=-77.9;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = -670.;    X_CUT_right = -610.;
    Y_CUT_left = -2048.;    Y_CUT_right = 2048.;
}

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

/////////** UA9 2015 - February **/////////

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

if( (runnumb >= 1720) &&
   (runnumb <= 1750) )
{
    deflXmar1=-300.;    deflXmar2=300.;
    deflYmar1=-100.;    deflYmar2=100.;
    
    thXin=64.;
    thXout=256.;
    
    torsion_posYthX=-9.65;//murad on mm
    torsion_posXthX=0.;//murad on mm
    meanchang = 80.;
    
    X_CUT_left = 170.;    X_CUT_right = 540.;
    Y_CUT_left = -4096.;    Y_CUT_right = 4096.;
}

/////////** UA9 2015 - May **/////////

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

if( (runnumb >= 2442) &&
   (runnumb <= 2442) )
{
    deflXmar1=-300.;    deflXmar2=300.;
    deflYmar1=-100.;    deflYmar2=100.;
    
    thXin=64.;
    thXout=256.;
    
    torsion_posYthX=10.;//25.;//murad on mm
    torsion_posXthX=0.;//murad on mm
    meanchang = -55.;
    
    X_CUT_left = -1430.;    X_CUT_right = 500.;
    Y_CUT_left = -4096.;    Y_CUT_right = 4096.;
}

if( (runnumb >= 2844) &&
   (runnumb <= 2850) )
{
    deflXmar1=-512.+100.;    deflXmar2=100.;
    deflYmar1=-64.;    deflYmar2=64.;
    
    thXin=32.;
    thXout=32.;
    
    torsion_posYthX=10.;//25.;//murad on mm
    torsion_posXthX=0.;//murad on mm
    meanchang = -55.;
    
    X_CUT_left = -200.;    X_CUT_right = 260.;
    Y_CUT_left = -4096.;    Y_CUT_right = 4096.;
}

if( (runnumb >= 2726) &&
   (runnumb <= 2732) )
{
    deflXmar1=-275;    deflXmar2=235.;
    deflYmar1=-128.;    deflYmar2=128.;
    
    thXin=64.;
    thXout=64.;
    
    torsion_posYthX=6.4+4.8;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -280.;
    
    X_CUT_left = -200.;    X_CUT_right = 300.;
    Y_CUT_left = -4096.;    Y_CUT_right = 4096.;
}

if( (runnumb >= 2753) &&
   (runnumb <= 2756) )
{
    deflXmar1=-128.;    deflXmar2=128.;
    deflYmar1=-128.;    deflYmar2=128.;
    
    thXin=32.;
    thXout=32.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -50.;
    
    X_CUT_left = -200.;    X_CUT_right = 300.;
    Y_CUT_left = 1000.;    Y_CUT_right = 4096.;
}

if( (runnumb >= 3078) &&
   (runnumb <= 3080))
{
    deflXmar1=-128.;    deflXmar2=128.;
    deflYmar1=-128.;    deflYmar2=128.;
    
    thXin=32.;
    thXout=32.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -50.;
    
    X_CUT_left = -500.;    X_CUT_right = 1000.;
    Y_CUT_left = 1000.;    Y_CUT_right = 4096.;
}

/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////
/////////*************************/////////

//PHOTAG/CHANEL 2015


if( (runnumb >= 110000) &&
   (runnumb <= 120000))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -168.;
    
    X_CUT_left = 7700. - 415.;    X_CUT_right = 7700. + 415.;
    Y_CUT_left = 11000.-500.;    Y_CUT_right = 11000.+500.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}


if( (runnumb >= 110000) &&
   (runnumb <= 110129))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=256.;
    thXout=512.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = 0;    X_CUT_right = 18000;
    Y_CUT_left = 0;    Y_CUT_right = 18000.;
}

if( (runnumb >= 110130) &&
   (runnumb <= 110132))
{
    deflXmar1=-256.;    deflXmar2=128.;
    deflYmar1=-128.;    deflYmar2=128.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = 9900;    X_CUT_right = 10700;
    Y_CUT_left = 4000;    Y_CUT_right = 14000.;
}

if( (runnumb >= 110133) &&
   (runnumb <= 110138))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = 10032.-450.;    X_CUT_right = 10032.+450.;
    Y_CUT_left = 9000.-1000.;    Y_CUT_right = 9000.+1000.;
}

if( (runnumb >= 110139) &&
   (runnumb <= 110160))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    
    X_CUT_left = 0;    X_CUT_right = 18000;
    Y_CUT_left = 0;    Y_CUT_right = 18000;
}
if( (runnumb >= 110161) &&
   (runnumb <= 110178))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -168.;
    
    X_CUT_left = 7700. - 400.;    X_CUT_right = 7700. + 400.;
    Y_CUT_left = 11000.-200.;    Y_CUT_right = 11000.+200.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}



if( (runnumb >= 112200) &&
   (runnumb <= 112328))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -168.;
    
    X_CUT_left = 11635.5 + 77.5 - 690.;    X_CUT_right = 11635.5 + 77.5 + 570.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

if( (runnumb >= 112329) &&
   (runnumb <= 112353))
{
    deflXmar1=-2048.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=256.;
    thXout=256.;
    
    torsion_posYthX=-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -168.;
    
    X_CUT_left = 11163. - 150.;    X_CUT_right = 11163. +150.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

if( (runnumb >= 112354) &&
   (runnumb <= 112432))
{
    deflXmar1=-256.;    deflXmar2=256.;
    deflYmar1=-256.;    deflYmar2=256.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    
    X_CUT_left = 11635.5 + 77.5 - 900.;    X_CUT_right = 11635.5 + 77.5 + 800.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

if( (runnumb >= 112432) &&
   (runnumb <= 112454))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    
    X_CUT_left = 11635.5 + 77.5 -700 -520. - 900.;    X_CUT_right = 11635.5 + 77.5 -700 -520. + 800.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

if( (runnumb >= 112455 && runnumb <= 112456 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    
    X_CUT_left = 11635.5 + 77.5 -700 -520. - 900.;    X_CUT_right = 11635.5 + 77.5 -700 -520. + 800.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

if( (runnumb >= 112457 && runnumb <= 112487 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=-13.; //-40.+6.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 11635.5 + 77.5 -700 -520. - 900.;    X_CUT_right = 11635.5 + 77.5 -700 -520. + 800.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}
if( (runnumb >= 112488 && runnumb <= 112517 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=-13.;//-40.+6.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 9400.;    X_CUT_right = 11000.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}
//BACK to positron Ge <110> axis (100) plane
if( (runnumb >= 112530 && runnumb <= 112554 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=0.;//-13.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = -100.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 9750.;    X_CUT_right = 11250.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

//PWO 23 Luglio positrons
if( (runnumb >= 112555 && runnumb <= 112582 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 9750.;    X_CUT_right = 11250.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

//PWO 25 Luglio electrons
if( (runnumb >= 112583 && runnumb <= 112592  ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 9250.;    X_CUT_right = 11040.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}

//Multi-strip sabbiata 26 Luglio positrons
if( (runnumb >= 112594 && runnumb <= 112600 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 10290.;    X_CUT_right = 10690.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}
//ST109 Si-<110>- new holder 26 Luglio positrons
if( (runnumb >= 112601 && runnumb <= 113000 ))
{
    deflXmar1=-512.;    deflXmar2=512.;
    deflYmar1=-512.;    deflYmar2=512.;
    
    thXin=512.;
    thXout=256.;
    
    torsion_posYthX=0.;//murad on mm
    torsion_posXthX=0.;//20.;//murad on mm
    meanchang = 0.;
    //X_CUT_left = 10700.;    X_CUT_right = 11250;
    X_CUT_left = 10040.;    X_CUT_right = 10750.;
    Y_CUT_left = 10000. - 4000.;    Y_CUT_right = 10000. + 4000.;
    //Y_CUT_left = 10500;    Y_CUT_right = 11500;
    
    DEVA_CUT = 600.;
}
