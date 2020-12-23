function [D, SL, out] = AssignDS_v2(QD_pdf,DS_pdf,Qr)
dmin = max([min(QD_pdf.Y), min(DS_pdf.X), 1.699]);
dmax = min(max(QD_pdf.Y), max(DS_pdf.X));
slmin = max(min(DS_pdf.Y), 1.699);
slmax = min(max(DS_pdf.Y), 3.1);

out = false;
cnt = 0;
while true
   cnt = cnt + 1;
   D = dmin + (dmax - dmin)*rand;
   SL = slmin + (slmax - slmin)*rand;
   r1 = rand;
   r2 = rand;
   R1 = QD_pdf.F(Qr, D);
   R2 = DS_pdf.F(D, SL);
   if r1 < R1 && r2 < R2
       out = true;
       break;
   end
   if cnt > 500
       break;
   end
end


