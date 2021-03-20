\\default(parisizemax,100m);
g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;
nval = 682492462409094395392022581537473179285250139967739310024802121913471471;


\\La stratégie est simple : Pohlig Hellman en ayant recours à BSGS

\\========================================================
\\                Baby Step Giant Step
\\========================================================

babyshark(n,gen,val) = {
	m = ceil(sqrt(n));
	a = Map();
	mul = Mod(1,gen.mod);
	for(i=0,m,mapput(a,lift(mul),i);mul = mul*gen;);
	b = val;
	inv = gen^(-m);
	for(i=0, m, if(mapisdefined(a,lift(b),&j),return(i*m +j));b *= inv;);
}

\\print(babyshark(nval,g,A));
\\il se trouve que la valeur de nval est trop élevé pour faire directement un BSGS dessus. Nous allons donc avoir recours à l'algorithme de Pohlig-Hellman

\\=======================================================
\\		  Pohlig-Hellman
\\=======================================================
smallPH(geni,Ai,f1i,f2i) ={
\\inspiré du PH_prime_order de Julien. Mon code continuait à avoir des overflow du fait du manque de RAM
	x=0;
	nwgen = geni^(f1i^(f2i-1));
	for(j=0, f2i-1, Aj = (geni^(-x)*Ai)^(f1i^(f2i-1-j));d= babyshark(f1i,nwgen,Aj);x+=f1i^j*d;);
	return(x);
}
phi = eulerphi(nval);
f = factor(phi);
a2 = vector(4);
for(i=1,4,nvali = f[i,1]^f[i,2];geni =g^(phi/nvali);Ai = Mod(A,nval)^(phi/nvali);a2[i]=Mod(smallPH(geni,Ai,f[i,1],f[i,2]),nvali););
print(lift(chinese(a2)));

