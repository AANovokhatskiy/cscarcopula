#include "generated_pdf.h" 
 
std::vector<double> pdf(const std::vector<std::vector<double>>& u, const double r)
{
    size_t n = u.size();
    std::vector<double> result(n);
 
    for (size_t i = 0; i < n; ++i) {
        double u0 = u[i][0];
        double u1 = u[i][1];
        double u2 = u[i][2];
        double u3 = u[i][3];
        double u4 = u[i][4];
        double u5 = u[i][5];
        double u6 = u[i][6];
        double u7 = u[i][7];
        double u8 = u[i][8];
        double tmp0 = log(u0);
        double tmp1 = log(u1);
        double tmp2 = log(u2);
        double tmp3 = log(u3);
        double tmp4 = log(u4);
        double tmp5 = log(u5);
        double tmp6 = log(u6);
        double tmp7 = log(u7);
        double tmp8 = log(u8);
        double tmp9 = pow(-tmp0, r);
        double tmp10 = pow(-tmp1, r);
        double tmp11 = pow(-tmp2, r);
        double tmp12 = pow(-tmp3, r);
        double tmp13 = pow(-tmp4, r);
        double tmp14 = pow(-tmp5, r);
        double tmp15 = pow(-tmp6, r);
        double tmp16 = pow(-tmp7, r);
        double tmp17 = pow(-tmp8, r);
        double tmp18 = tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
        double tmp19 = 1.0/r;
        double tmp20 = pow(tmp18, tmp19);
        double tmp21 = pow(r, 2);
        double tmp22 = 546*tmp21;
        double tmp23 = pow(r, 3);
        double tmp24 = pow(r, 4);
        double tmp25 = 22449*tmp24;
        double tmp26 = pow(r, 5);
        double tmp27 = pow(r, 6);
        double tmp28 = 118124*tmp27;
        double tmp29 = pow(r, 7);
        double tmp30 = pow(tmp18, 2*tmp19);
        double tmp31 = pow(tmp18, 3*tmp19);
        double tmp32 = pow(tmp18, 4*tmp19);
        double tmp33 = pow(tmp18, 5*tmp19);
        double tmp34 = pow(tmp18, 6*tmp19);
        double tmp35 = pow(tmp18, 7*tmp19);
        result[i] = tmp10*tmp11*tmp12*tmp13*tmp14*tmp15*tmp16*tmp17*pow(tmp18, tmp19 - 9)*tmp9*(-40320*pow(r, 8) - 4572*r*tmp20 + 34776*r*tmp30 - 61236*r*tmp31 + 37800*r*tmp32 - 9576*r*tmp33 + 1008*r*tmp34 - 36*r*tmp35 + 36*r - pow(tmp18, 8*tmp19) + 34398*tmp20*tmp21 - 140616*tmp20*tmp23 + 336735*tmp20*tmp24 - 470988*tmp20*tmp26 + 354372*tmp20*tmp27 - 109584*tmp20*tmp29 + 255*tmp20 - 164346*tmp21*tmp30 + 191100*tmp21*tmp31 - 76440*tmp21*tmp32 + 11466*tmp21*tmp33 - tmp22*tmp34 - tmp22 + 408240*tmp23*tmp30 - 294840*tmp23*tmp31 + 68040*tmp23*tmp32 - 4536*tmp23*tmp33 + 4536*tmp23 - 561225*tmp24*tmp30 + 224490*tmp24*tmp31 - tmp25*tmp32 - tmp25 + 403704*tmp26*tmp30 - 67284*tmp26*tmp31 + 67284*tmp26 - tmp28*tmp30 - tmp28 + 109584*tmp29 - 3025*tmp30 + 7770*tmp31 - 6951*tmp32 + 2646*tmp33 - 462*tmp34 + 36*tmp35 - 1)*exp(-tmp20)/(tmp0*tmp1*tmp2*tmp3*tmp4*tmp5*tmp6*tmp7*tmp8*u0*u1*u2*u3*u4*u5*u6*u7*u8);
    }
  return result;
}

std::vector<double> pdf(const std::vector<std::vector<double>>& u, const std::vector<double>& r)
{
    size_t n = u.size();
    std::vector<double> result(n);
 
    for (size_t i = 0; i < n; ++i) {
        double u0 = u[i][0];
        double u1 = u[i][1];
        double u2 = u[i][2];
        double u3 = u[i][3];
        double u4 = u[i][4];
        double u5 = u[i][5];
        double u6 = u[i][6];
        double u7 = u[i][7];
        double u8 = u[i][8];
        double tmp0 = log(u0);
        double tmp1 = log(u1);
        double tmp2 = log(u2);
        double tmp3 = log(u3);
        double tmp4 = log(u4);
        double tmp5 = log(u5);
        double tmp6 = log(u6);
        double tmp7 = log(u7);
        double tmp8 = log(u8);
        double tmp9 = pow(-tmp0, r[i]);
        double tmp10 = pow(-tmp1, r[i]);
        double tmp11 = pow(-tmp2, r[i]);
        double tmp12 = pow(-tmp3, r[i]);
        double tmp13 = pow(-tmp4, r[i]);
        double tmp14 = pow(-tmp5, r[i]);
        double tmp15 = pow(-tmp6, r[i]);
        double tmp16 = pow(-tmp7, r[i]);
        double tmp17 = pow(-tmp8, r[i]);
        double tmp18 = tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
        double tmp19 = 1.0/r[i];
        double tmp20 = pow(tmp18, tmp19);
        double tmp21 = pow(r[i], 2);
        double tmp22 = 546*tmp21;
        double tmp23 = pow(r[i], 3);
        double tmp24 = pow(r[i], 4);
        double tmp25 = 22449*tmp24;
        double tmp26 = pow(r[i], 5);
        double tmp27 = pow(r[i], 6);
        double tmp28 = 118124*tmp27;
        double tmp29 = pow(r[i], 7);
        double tmp30 = pow(tmp18, 2*tmp19);
        double tmp31 = pow(tmp18, 3*tmp19);
        double tmp32 = pow(tmp18, 4*tmp19);
        double tmp33 = pow(tmp18, 5*tmp19);
        double tmp34 = pow(tmp18, 6*tmp19);
        double tmp35 = pow(tmp18, 7*tmp19);
        result[i] = tmp10*tmp11*tmp12*tmp13*tmp14*tmp15*tmp16*tmp17*pow(tmp18, tmp19 - 9)*tmp9*(-40320*pow(r[i], 8) - 4572*r[i]*tmp20 + 34776*r[i]*tmp30 - 61236*r[i]*tmp31 + 37800*r[i]*tmp32 - 9576*r[i]*tmp33 + 1008*r[i]*tmp34 - 36*r[i]*tmp35 + 36*r[i] - pow(tmp18, 8*tmp19) + 34398*tmp20*tmp21 - 140616*tmp20*tmp23 + 336735*tmp20*tmp24 - 470988*tmp20*tmp26 + 354372*tmp20*tmp27 - 109584*tmp20*tmp29 + 255*tmp20 - 164346*tmp21*tmp30 + 191100*tmp21*tmp31 - 76440*tmp21*tmp32 + 11466*tmp21*tmp33 - tmp22*tmp34 - tmp22 + 408240*tmp23*tmp30 - 294840*tmp23*tmp31 + 68040*tmp23*tmp32 - 4536*tmp23*tmp33 + 4536*tmp23 - 561225*tmp24*tmp30 + 224490*tmp24*tmp31 - tmp25*tmp32 - tmp25 + 403704*tmp26*tmp30 - 67284*tmp26*tmp31 + 67284*tmp26 - tmp28*tmp30 - tmp28 + 109584*tmp29 - 3025*tmp30 + 7770*tmp31 - 6951*tmp32 + 2646*tmp33 - 462*tmp34 + 36*tmp35 - 1)*exp(-tmp20)/(tmp0*tmp1*tmp2*tmp3*tmp4*tmp5*tmp6*tmp7*tmp8*u0*u1*u2*u3*u4*u5*u6*u7*u8);
    }
  return result;
}

double pdf(const std::vector<double>& u, const double r)
{
    double result = 0.0;
 
    double u0 = u[0];
    double u1 = u[1];
    double u2 = u[2];
    double u3 = u[3];
    double u4 = u[4];
    double u5 = u[5];
    double u6 = u[6];
    double u7 = u[7];
    double u8 = u[8];
    double tmp0 = log(u0);
    double tmp1 = log(u1);
    double tmp2 = log(u2);
    double tmp3 = log(u3);
    double tmp4 = log(u4);
    double tmp5 = log(u5);
    double tmp6 = log(u6);
    double tmp7 = log(u7);
    double tmp8 = log(u8);
    double tmp9 = pow(-tmp0, r);
    double tmp10 = pow(-tmp1, r);
    double tmp11 = pow(-tmp2, r);
    double tmp12 = pow(-tmp3, r);
    double tmp13 = pow(-tmp4, r);
    double tmp14 = pow(-tmp5, r);
    double tmp15 = pow(-tmp6, r);
    double tmp16 = pow(-tmp7, r);
    double tmp17 = pow(-tmp8, r);
    double tmp18 = tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
    double tmp19 = 1.0/r;
    double tmp20 = pow(tmp18, tmp19);
    double tmp21 = pow(r, 2);
    double tmp22 = 546*tmp21;
    double tmp23 = pow(r, 3);
    double tmp24 = pow(r, 4);
    double tmp25 = 22449*tmp24;
    double tmp26 = pow(r, 5);
    double tmp27 = pow(r, 6);
    double tmp28 = 118124*tmp27;
    double tmp29 = pow(r, 7);
    double tmp30 = pow(tmp18, 2*tmp19);
    double tmp31 = pow(tmp18, 3*tmp19);
    double tmp32 = pow(tmp18, 4*tmp19);
    double tmp33 = pow(tmp18, 5*tmp19);
    double tmp34 = pow(tmp18, 6*tmp19);
    double tmp35 = pow(tmp18, 7*tmp19);
        result += tmp10*tmp11*tmp12*tmp13*tmp14*tmp15*tmp16*tmp17*pow(tmp18, tmp19 - 9)*tmp9*(-40320*pow(r, 8) - 4572*r*tmp20 + 34776*r*tmp30 - 61236*r*tmp31 + 37800*r*tmp32 - 9576*r*tmp33 + 1008*r*tmp34 - 36*r*tmp35 + 36*r - pow(tmp18, 8*tmp19) + 34398*tmp20*tmp21 - 140616*tmp20*tmp23 + 336735*tmp20*tmp24 - 470988*tmp20*tmp26 + 354372*tmp20*tmp27 - 109584*tmp20*tmp29 + 255*tmp20 - 164346*tmp21*tmp30 + 191100*tmp21*tmp31 - 76440*tmp21*tmp32 + 11466*tmp21*tmp33 - tmp22*tmp34 - tmp22 + 408240*tmp23*tmp30 - 294840*tmp23*tmp31 + 68040*tmp23*tmp32 - 4536*tmp23*tmp33 + 4536*tmp23 - 561225*tmp24*tmp30 + 224490*tmp24*tmp31 - tmp25*tmp32 - tmp25 + 403704*tmp26*tmp30 - 67284*tmp26*tmp31 + 67284*tmp26 - tmp28*tmp30 - tmp28 + 109584*tmp29 - 3025*tmp30 + 7770*tmp31 - 6951*tmp32 + 2646*tmp33 - 462*tmp34 + 36*tmp35 - 1)*exp(-tmp20)/(tmp0*tmp1*tmp2*tmp3*tmp4*tmp5*tmp6*tmp7*tmp8*u0*u1*u2*u3*u4*u5*u6*u7*u8);
  return result;
}

double transform(const double r)
{
    double result = 0.0;
  result = fmin(40, pow(r, 2) + 1);
  return result;
}

std::vector<double> transform(const std::vector<double>& r)
{
    size_t n = r.size();
    std::vector<double> result(n);
 
    for (size_t i = 0; i < n; ++i) {
        result[i] = fmin(40, pow(r[i], 2) + 1);
    }
  return result;
}
