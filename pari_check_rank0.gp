\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\ pari_check_rank0.gp -- Rank-0 certification for E_A and E'_A
\\\
\\\ PURPOSE
\\\ -------
\\\ This script uses PARI/GP's ellrank() (Simon two-descent) to certify
\\\ that certain specialisations of the elliptic curves E_A and E'_A have
\\\ Mordell--Weil rank 0 over Q.
\\\
\\\ For each coprime pair (a, b) with 1 <= a < b <= 19 we set the
\\\ rational parameter s = a/b and compute:
\\\
\\\   E_A  :  y^2 = x^3 + A*x^2 - 4*x - 4*A
\\\            where A = 2 - 4*c^2,  c = (s^4 - 6s^2 + 1)/(1+s^2)^2.
\\\
\\\   E'_A :  y^2 = x^3 + (2/V)*x^2 - x - (2/V)
\\\            where V = 8s(1-s^2)/(1+s^2)^2.
\\\
\\\ ellrank() returns a pair [lower_bound, upper_bound] for the rank.
\\\ An upper bound of 0 certifies rank = 0; no further Mordell--Weil
\\\ search is needed.
\\\
\\\ PAPER CLAIM SUPPORTED
\\\ ---------------------
\\\ "PARI/GP's ellrank certifies rk E_A = 0 for 42 values and
\\\  rk E'_A = 0 for 54 values (upper bound 0 from Simon 2-descent)."
\\\
\\\ The script prints a summary showing the number of certified rank-0
\\\ specialisations for each family.  The expected output is:
\\\   E_A  : 42 certified rank-0
\\\   E'_A : 54 certified rank-0
\\\
\\\ USAGE
\\\ -----
\\\   gp < pari_check_rank0.gp
\\\   (requires PARI/GP with the ellrank / ell.gp package)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ea_cert = 0; ea_total = 0; ea_rk1 = 0; ea_unsicher = 0;
ep_cert = 0; ep_total = 0; ep_rk1 = 0; ep_unsicher = 0;

print("=== E_A: y^2 = x^3 + A*x^2 - 4*x - 4*A ===");
print("");

for(a_val = 1, 19,
  for(b_val = a_val+1, 19,
    if(gcd(a_val, b_val) == 1,
      s0 = a_val/b_val;
      c0 = (s0^4 - 6*s0^2 + 1)/(1+s0^2)^2;
      A0 = 2 - 4*c0^2;

      E = ellinit([0, A0, 0, -4, -4*A0], flag=1);
      ea_total++;

      rk = ellrank(E);
      lb = rk[1]; ub = rk[2];

      if(ub == 0,
        ea_cert++;
        if(ea_cert <= 10,
          print("  CERT r=0: s=", a_val, "/", b_val, "  rank=[", lb, ",", ub, "]")
        )
      ,
        if(lb > 0,
          ea_rk1++;
        ,
          ea_unsicher++;
          if(ea_unsicher <= 5,
            print("  UNSICHER: s=", a_val, "/", b_val, "  rank=[", lb, ",", ub, "]")
          )
        )
      )
    )
  )
);

print("");
print("E_A: Getestet=", ea_total, "  Zertifiziert_r0=", ea_cert, "  Rang>0=", ea_rk1, "  Unsicher=", ea_unsicher);

print("");
print("=== E'_A ===");
print("");

for(a_val = 1, 19,
  for(b_val = a_val+1, 19,
    if(gcd(a_val, b_val) == 1,
      s0 = a_val/b_val;
      V0 = 8*s0*(1-s0^2)/(1+s0^2)^2;
      if(V0 != 0,
        coeff = 2/V0;

        E = ellinit([0, coeff, 0, -1, -coeff], flag=1);
        ep_total++;

        rk = ellrank(E);
        lb = rk[1]; ub = rk[2];

        if(ub == 0,
          ep_cert++;
          if(ep_cert <= 10,
            print("  CERT r=0: s=", a_val, "/", b_val, "  rank=[", lb, ",", ub, "]")
          )
        ,
          if(lb > 0,
            ep_rk1++;
          ,
            ep_unsicher++;
            if(ep_unsicher <= 5,
              print("  UNSICHER: s=", a_val, "/", b_val, "  rank=[", lb, ",", ub, "]")
            )
          )
        )
      )
    )
  )
);

print("");
print("E'_A: Getestet=", ep_total, "  Zertifiziert_r0=", ep_cert, "  Rang>0=", ep_rk1, "  Unsicher=", ep_unsicher);

print("");
print("=== FAZIT ===");
if(ea_cert > 0, print("E_A:  OK (", ea_cert, " zertifiziert)"), print("E_A:  PROBLEM!"));
if(ep_cert > 0, print("E'_A: OK (", ep_cert, " zertifiziert)"), print("E'_A: PROBLEM!"));
}
quit;
