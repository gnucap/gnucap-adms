// these devices are provided by the simulator and need
// not be implemented

extern module R(p,n);
   inout p,n;
	electrical p,n;

	parameter real r=1e3 from [0:inf);
endmodule

extern module C(p,n);
   inout p,n;
	electrical p,n;

	parameter real c=1e-6 from [0:inf);
endmodule

extern module rcd_exp(p, n, e);
	inout p, n;
	output e;
	electrical p, n;
	degradational e;

	parameter real weight=1 from (0:inf);
	parameter real re1=1e3 from [0:inf);
	parameter real re0=1e3 from (-inf:inf);
	parameter real rc1=1e3 from (-inf:0];
	parameter real rc0=1e3 from (-inf:inf);
endmodule
