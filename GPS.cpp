#define C_EARTH (double)6378137.0
#define C_PI (double)3.141592653589793
#define DEG2RAD(DEG) ((DEG) * ((C_PI) / (180.0)))

typedef struct
{
    double lat;
    double lon;
} NavSatFix;

//已知两点gps_position和start_gps_position,求两点间直线距离
double local_Offset(NavSatFix& gps_position, NavSatFix& start_gps_position)
{
  double radLat1 = DEG2RAD(gps_position.lat);
  double radLat2 = DEG2RAD(start_gps_position.lat);
  double a = radLat1 - radLat2;
  double b = DEG2RAD(gps_position.lon) - DEG2RAD(start_gps_position.lon);
  double dst = 2 * asin((sqrt(pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2))));
  dst = dst * C_EARTH;
  return round(dst * 10000) / 10000;
}

//已知一点start_gps_position,距离d(单位m), 角度ang(北0度，东90度，南180度，西270度)，求另一点经纬度
NavSatFix fun(NavSatFix& start_gps_position, double d, double ang)
{
    NavSatFix data;
    data.lat = start_gps_position.lat + d * cos(ang * C_PI /180.0) / 111122.1976989967;
    data.lon = start_gps_position.lon + d * sin(ang * C_PI /180.0) / (111122.1976989967 * cos(start_gps_position.lat * C_PI /180));
    return data;
}
