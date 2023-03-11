// note when we do spherical coordinates (r,theta, phi) that r and theta are just like polar coordinates on the XY plane 
// and phi is the angle from the positive Z axis. All angles are in radians.
// for 38 degrees try 25 and 32.
var sphRadList = []; // (r,theta, phi) in radians;

function lower() {
  document.getElementById("Urise").value = document.getElementById("Lrise").value;
  document.getElementById("Urun").value = document.getElementById("Lrun").value;
  init();
}

function miser() {
  if (document.getElementById("miser").checked) {
    document.getElementById("upper").style.display = "none";
    document.getElementById("note").style.display = "block";

  } else {
    document.getElementById("upper").style.display = "block";
    document.getElementById("Urise").value=document.getElementById("Lrise").value;
    document.getElementById("Urun").value=document.getElementById("Lrun").value;
    document.getElementById("note").style.display = "none";
  }
  init();
}

// this reads changes user made to angle & elevation textarea.
function userInput() {
  var myTextarea = "["+document.getElementById("angElev").value+"]";
  var angElev = JSON.parse(myTextarea);
  sphRadList = [];
  angElev.forEach(function(pt) {
    var theta = deg2rad(pt[0]);
    var phi = deg2rad(90-pt[1]);
    sphRadList.push([1,theta,phi]);
  });
  if (document.getElementById("miser").checked) {
    var lastPt = angElev[angElev.length-1];
    if (lastPt[0] != 0) {
      var phi = deg2rad(90-lastPt[1]);
      sphRadList.push([1,0,phi]);
    }
    for (var i = angElev.length-2; i > -1; i--) { 
      var theta = -deg2rad(angElev[i][0]);
      var phi = deg2rad(90-angElev[i][1]);
      sphRadList.push([1,theta,phi]);
    }
  }
//alert(sphRadList);
  initB();
}

function deg2rad(deg) {
  return(deg*Math.PI/180);
}

function rad2deg(rad) {
  return (rad*180/Math.PI);
}

function sph2cart(sph) {
  var r = sph[0];
  var theta = sph[1];
  var phi = sph[2];
  var xx = r*Math.sin(phi)*Math.cos(theta);
  var yy = r*Math.sin(phi)*Math.sin(theta);
  var zz = r*Math.cos(phi);
  return([xx,yy,zz]);
}

function cart2sph(cart) {
  var x = cart[0];
  var y = cart[1];
  var z = cart[2];
  var theta, phi;
  var r = Math.sqrt(x*x+y*y+z*z);
  if (r===0) {
    theta = 0;
    phi = 0;
  } else {
    phi = Math.acos(z/r);
    if (x===0) {
      if (y>0) {
        theta = Math.PI/2;
      } else if (y<0) {
        theta = 3*Math.PI/2;
      } else {
        theta = 0;
      }
    } else { // if x != 0
      theta = Math.atan(y/x);
      if (x<0) {theta += Math.PI;}
    }
  } // end r != 0
  return([r,theta,phi]);
} // end cart2sph()

function dot(vect1, vect2) {
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}

function cross(vect1,vect2) {
  var newX = vect1[1]*vect2[2]-vect1[2]*vect2[1];
  var newY = vect1[2]*vect2[0]-vect1[0]*vect2[2];
  var newZ = vect1[0]*vect2[1]-vect1[1]*vect2[0];
  return([newX,newY,newZ]);
}

function vectDiff(vect1, vect2) {
  return([vect1[0]-vect2[0],vect1[1]-vect2[1],vect1[2]-vect2[2]]);
}

function normalize(vect) {
  var vectLen = Math.sqrt(dot(vect,vect));
  if (vectLen===0) {return ([0,0,0]);}
  else {return([vect[0]/vectLen,vect[1]/vectLen,vect[2]/vectLen]);}
}

function uVectAngle(uVect1,uVect2) {
  return(Math.acos(dot(uVect1,uVect2)));
}

function rotMatrix(axisVect,ang) {
  var xx = axisVect[0];
  var yy = axisVect[1];
  var zz = axisVect[2];
  var r11=Math.cos(ang) + xx*xx*(1-Math.cos(ang));
  var r21=yy*xx*(1-Math.cos(ang)) + zz*Math.sin(ang);
  var r31=zz*xx*(1-Math.cos(ang)) - yy*Math.sin(ang);
  var r12=xx*yy*(1-Math.cos(ang)) - zz*Math.sin(ang);
  var r22=Math.cos(ang) + yy*yy*(1-Math.cos(ang));
  var r32=zz*yy*(1-Math.cos(ang)) + xx*Math.sin(ang);
  var r13=xx*zz*(1-Math.cos(ang)) + yy*Math.sin(ang);
  var r23=yy*zz*(1-Math.cos(ang)) - xx*Math.sin(ang);
  var r33=Math.cos(ang) + zz*zz*(1-Math.cos(ang));
  return([[r11,r21,r31],[r12,r22,r32],[r13,r23,r33]]);
}

function multMatVect(matrix,vect) {
  var newX = vect[0]*matrix[0][0]+vect[1]*matrix[0][1]+vect[2]*matrix[0][2];
  var newY = vect[0]*matrix[1][0]+vect[1]*matrix[1][1]+vect[2]*matrix[1][2];
  var newZ = vect[0]*matrix[2][0]+vect[1]*matrix[2][1]+vect[2]*matrix[2][2];
  return([newX,newY,newZ]);
}

function findArea(mySph) {
  var myCart = [];
  var triangles = [];
  var totalArea = 0;
  mySph.forEach(function(step) {
    var cart = sph2cart([1,step[1],step[2]]);
    myCart.push(cart);
  });

  for (var i=2;i<myCart.length;i++) {
    triangles.push([myCart[0],myCart[i-1],myCart[i]]);
  }
//alert(JSON.stringify(triangles));
  triangles.forEach(function(tri) {
    var side1 = uVectAngle(tri[0],tri[1]);
    var side2 = uVectAngle(tri[0],tri[2]);
    var side3 = uVectAngle(tri[1],tri[2]);
    var s = (side1+side2+side3)/2;
    var tans = Math.tan(.5*s)*Math.tan(.5*(s-side1))*Math.tan(.5*(s-side2))*Math.tan(.5*(s-side3));
    if (tans<0) {tans=0};
    var area= 4*Math.atan(Math.sqrt(tans));
    if (area>0.000001) {
      var middle = [tri[0][0]+tri[1][0]+tri[2][0],tri[0][1]+tri[1][1]+tri[2][1],tri[0][2]+tri[1][2]+tri[2][2]];
      middle = normalize(middle);
      var CCW=dot(cross(vectDiff(tri[1],tri[0]),vectDiff(tri[2],tri[1])),tri[1]);
      if (CCW<0) {area=-area}
    }
    totalArea += area;
  });
  return(totalArea);
}

function init() {
  if (document.getElementById("miser").checked) {
    initMiser();
    initB();
  } else {
    initNoMiser();
    initB();
  }
}

function initB() {
  // set angle & elevation 
  var angElev2 = "";
  if (document.getElementById("miser").checked) {
    for (var i=0;i<sphRadList.length/2;i++) {
      angElev2 = angElev2.concat("[" 
         + Math.round(100*rad2deg(sphRadList[i][1]))/100+"," 
         + Math.round(100*(90-rad2deg(sphRadList[i][2])))/100+"],\r\n");
    }
  } else {
    sphRadList.forEach(function(pt) {
      angElev2 = angElev2.concat("[" 
         + Math.round(100*rad2deg(pt[1]))/100+"," 
         + Math.round(100*(90-rad2deg(pt[2])))/100+"],\r\n");
    });
  }
  angElev2 = angElev2.substr(0,angElev2.length-3);
  document.getElementById("angElev").value = angElev2;

  // set area and area2
  var standSphRadList = [];
  standSphRadList.push(sphRadList[0]);
  standSphRadList.push([1,sphRadList[0][1],Math.PI/2]);
  standSphRadList.push([1,sphRadList[sphRadList.length-1][1],Math.PI/2]);
  standSphRadList.push(sphRadList[sphRadList.length-1]);
  document.getElementById("area").innerHTML = JSON.stringify(findArea(standSphRadList));
  document.getElementById("area2").innerHTML = JSON.stringify(findArea(sphRadList));

  //set cartesian vectors
  var cartList = []; //vector as cartesian (x,y,z);
  sphRadList.forEach(function(step) {
    var cart = sph2cart([1,step[1],step[2]]);
    cartList.push(cart);
  });

  // find length of sides of sphere polygon
  var mySides = "";
  for (var i=1;i<cartList.length;i++) {
    mySides = mySides.concat("" + Math.round(1000*uVectAngle(cartList[i-1],cartList[i]))/1000 +"\r\n");
  }
  document.getElementById("sides").innerHTML = mySides;

  // long process to find bevel and miter
  // find an orientation vector to start. Note that we aren't allowed to start with the z axis.
  var zAxis = [0,0,1];
  var olddir = cartList[0];
  var perpenDir = cross(zAxis,olddir);
  var oldori = normalize(cross(olddir,perpenDir));
  var cuts = [];
  for (var i=1;i<cartList.length;i++) {
    var dir = cartList[i];
    // oldper is perpendicular to olddir and oldori.;
    // (olddir, oldper, oldori) is a new coordinate system.;
    var oldper = cross(oldori,olddir);
    // n is dir in (olddir, oldper, oldori) coordinate system;
    var nx = dot(dir,olddir);
    var ny = dot(dir,oldper);
    var nz = dot(dir,oldori);
    var nSph = cart2sph([nx,ny,nz]);
    var phiDeg = rad2deg(nSph[2]);
    var thetaDeg = rad2deg(nSph[1]);

    // for bevel and miter, if you have the main stock on the left side of a miter saw, 
    // the bevel is positive if you vertically tilt the saw to the right like this: /.
    // the miter is positive if you horizontally angle the saw to the right like this: \.

    var bevel=Math.round((90-Math.abs(phiDeg))/2*100)/100;
    var miter=Math.round(thetaDeg/2*100)/100;
    if (phiDeg > 0) {bevel=-bevel;}
    // u is perpendicular to olddir and dir.;
    var u = cross(dir,olddir);
    u = normalize(u);
    // thet is the angle between olddir and dir;
    var thet = uVectAngle(dir,olddir);
    // rMatrix is the rotation matrix around vector u from olddir to dir;
    var rMatrix = rotMatrix(u,thet);
    // q should match dir. to check we did it right;
    var q = multMatVect(rMatrix, olddir);
    // this should be the new rotated orientation;
    var ori = multMatVect(rMatrix, oldori);
    cuts.push([Math.round(bevel*100)/100,Math.round(miter*100)/100]);
    oldori = ori;
    olddir = dir;
  } // end for loop

  var myCuts = "";
  if (document.getElementById("miser").checked) {
    for (var i = 1;i<cuts.length;i+=2) {
      cuts[i][0]*=-1;
      cuts[i][1]*=-1;
    }
    cuts.forEach(function(cut) {
      myCuts = myCuts.concat(" to " + cut +"\r\n");
      myCuts = myCuts.concat("" + (cut[0])+","+(cut[1]));
    });
  } else {
    cuts.forEach(function(cut) {
      myCuts = myCuts.concat(" to " + cut +"\r\n");
      myCuts = myCuts.concat("" + (-cut[0])+","+(-cut[1]));
    });
  }
  document.getElementById("bevMit").innerHTML = myCuts;

} //end initB()

// This will read off the upper and lower rise and run from the webpage and set the standard handrail
function initNoMiser() {
  sphRadList = [];

  var lowrise = parseFloat(document.getElementById("Lrise").value);
  var lowrun = parseFloat(document.getElementById("Lrun").value);
  var highrise = parseFloat(document.getElementById("Urise").value);
  var highrun = parseFloat(document.getElementById("Urun").value);
  // turnAngle is how far you turn to the left when going upstairs.
  var turnAngle = parseFloat(document.getElementById("Turn").value);

  //this is the lower stair;
  var theta = -deg2rad(turnAngle/2);
  var phi=Math.PI/2-Math.atan(lowrise/lowrun);
  sphRadList.push([1,theta,phi]);

  //this is the level part after the lower stair;
  phi=Math.PI/2;
  sphRadList.push([1,theta,phi]);

  //this is the level part before the upper stair;
  theta=-theta;
  sphRadList.push([1,theta,phi]);

  //this is the upper stair;
  phi=Math.PI/2-Math.atan(highrise/highrun);
  sphRadList.push([1,theta,phi]);
} // end initNoMiser


// This will read off the upper and lower rise and run from the webpage and set the standard miser handrail
function initMiser() {
  sphRadList = [];

  var lowrise = parseFloat(document.getElementById("Lrise").value);
  var lowrun = parseFloat(document.getElementById("Lrun").value);
  var highrise = lowrise;
  var highrun = lowrun;
  // turnAngle is how far you turn to the left when going upstairs.
  var turnAngle = parseFloat(document.getElementById("Turn").value);

  //this is the lower stair;
  var theta = -deg2rad(turnAngle/2);
  var phi=Math.PI/2-Math.atan(lowrise/lowrun);
  sphRadList.push([1,theta,phi]);

  //this is the level part after the lower stair;
  phi=Math.PI/2;
  sphRadList.push([1,theta,phi]);

  sphRadList.push([1,0,phi]);

  //this is the level part before the upper stair;
  theta=-theta;
  sphRadList.push([1,theta,phi]);

  //this is the upper stair;
  phi=Math.PI/2-Math.atan(highrise/highrun);
  sphRadList.push([1,theta,phi]);
} // end initMiser