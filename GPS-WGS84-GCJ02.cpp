
#define C_EARTH (double)6378137.0
#define A_EARTH (double)6378245.0
#define C_PI (double)3.141592653589793
#define EE (double)0.00669342162296594323
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



//判断该GCJ-02坐标点是否在国内，在国内则返回false    GCJ-02坐标为国内对WGS-84坐标加密、偏移后所得，仅国内使用GCJ-02坐标 高德地图为GCJ-02坐标
bool outOfChina(NavSatFix& gpsData) {
    if (gpsData.lon < 72.004 || gpsData.lon > 137.8347) return true;
    if (gpsData.lat < 0.8293 || gpsData.lat > 55.8271) return true;
    return false;
}

double transformLat(double x, double y) {
    double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * C_PI) + 20.0 * sin(2.0 * x * C_PI)) * 2.0 / 3.0;
    ret += (20.0 * sin(y * C_PI) + 40.0 * sin(y / 3.0 * C_PI)) * 2.0 / 3.0;
    ret += (160.0 * sin(y / 12.0 * C_PI) + 320 * sin(y * C_PI / 30.0)) * 2.0 / 3.0;
    return ret;
}

double transformLon(double x, double y) {
    double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * C_PI) + 20.0 * sin(2.0 * x * C_PI)) * 2.0 / 3.0;
    ret += (20.0 * sin(x * C_PI) + 40.0 * sin(x / 3.0 * C_PI)) * 2.0 / 3.0;
    ret += (150.0 * sin(x / 12.0 * C_PI) + 300.0 * sin(x / 30.0 * C_PI)) * 2.0 / 3.0;
    return ret;
}

//GCJ02坐标装换为WGS84坐标
NavSatFix GCJ02_To_WGS84(NavSatFix& gpsData) {
    if (outOfChina(gpsData)) return gpsData;
    double dLat = transformLat(gpsData.lon - 105.0, gpsData.lat - 35.0);
    double dLon = transformLon(gpsData.lon - 105.0, gpsData.lat - 35.0);
    double radLat = gpsData.lat / 180.0 * C_PI;
    double magic = sin(radLat);
    magic = 1 - EE * magic * magic;
    double sqrtMagic = sqrt(magic);
    dLat = (dLat * 180.0) / ((A_EARTH * (1 - EE)) / (magic * sqrtMagic) * C_PI);
    dLon = (dLon * 180.0) / (A_EARTH / sqrtMagic * cos(radLat) * C_PI);
    double mgLat = gpsData.lat + dLat;
    double mgLon = gpsData.lon + dLon;
    double latitude = gpsData.lat * 2 - mgLat;
    double lontitude = gpsData.lon * 2 - mgLon;
    return NavSatFix(latitude, lontitude);
}

//WGS84坐标装换为GCJ02坐标
NavSatFix WGS84_To_GCJ02(NavSatFix& gpsData) {
  if(outOfChina(gpsData)) {
    return gpsData;
  }
  double dLat = transformLat(gpsData.lon - 105.0, gpsData.lat - 35.0);
  double dLon = transformLon(gpsData.lon - 105.0, gpsData.lat - 35.0);
  double radLat = gpsData.lat / 180.0 * C_PI;
  double magic = sin(radLat);
  magic = 1 - EE * magic * magic;
  double sqrtMagic = sqrt(magic);
  dLat = (dLat * 180.0) / ((A_EARTH * (1 - EE)) / (magic * sqrtMagic) * C_PI);
  dLon = (dLon * 180.0) / (A_EARTH / sqrtMagic * cos(radLat) * C_PI);
  return NavSatFix(gpsData.lat + dLat, gpsData.lon + dLon);
}
