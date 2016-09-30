{ // contribution I(p,n)<+((V(p,n)/R));
	double e;
	double e_Vp_n;
	e_Vp_n = /*ep*/ ( /*pprobe*/ 1. /m->v_R); // 7672
	assert(is_number(e_Vp_n));
	double x0x8612a0 = e =  ((/* treval */  BP(p, n))/m->v_R); // ...
		assert(is_number(x0x8612a0)); /*tree7704*/

	{
		double vx = STATIC_RHS * e;
		assert(is_number(vx));
			_br_p_n_xxx[0] += vx;
			// contribute(p, n, 0, vx);
	}

	{ // ?<+ V in V(p,n) no
		double e = e_Vp_n;
		double ceq;
		ceq = - e * BR(p, n);
		contribute(_br_p_n_xxx, 0,  ceq); // offset
		contribute(_br_p_n_xxx, P__br_p_n__br_p_n_V, e);
	}

} // contribution
