/* mpfr_gmp -- Limited gmp-impl emulator
   Modified version of the GMP files.

Copyright 2004-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

#ifndef MPFR_HAVE_GMP_IMPL

const struct bases mpfr_bases[257] =
{
  /*  0 */ {0.0},
  /*  1 */ {1e37},
  /*  2 */ {1.0000000000000000},
  /*  3 */ {0.6309297535714574},
  /*  4 */ {0.5000000000000000},
  /*  5 */ {0.4306765580733931},
  /*  6 */ {0.3868528072345416},
  /*  7 */ {0.3562071871080222},
  /*  8 */ {0.3333333333333333},
  /*  9 */ {0.3154648767857287},
  /* 10 */ {0.3010299956639812},
  /* 11 */ {0.2890648263178878},
  /* 12 */ {0.2789429456511298},
  /* 13 */ {0.2702381544273197},
  /* 14 */ {0.2626495350371935},
  /* 15 */ {0.2559580248098155},
  /* 16 */ {0.2500000000000000},
  /* 17 */ {0.2446505421182260},
  /* 18 */ {0.2398124665681314},
  /* 19 */ {0.2354089133666382},
  /* 20 */ {0.2313782131597592},
  /* 21 */ {0.2276702486969530},
  /* 22 */ {0.2242438242175754},
  /* 23 */ {0.2210647294575037},
  /* 24 */ {0.2181042919855316},
  /* 25 */ {0.2153382790366965},
  /* 26 */ {0.2127460535533632},
  /* 27 */ {0.2103099178571525},
  /* 28 */ {0.2080145976765095},
  /* 29 */ {0.2058468324604344},
  /* 30 */ {0.2037950470905062},
  /* 31 */ {0.2018490865820999},
  /* 32 */ {0.2000000000000000},
  /* 33 */ {0.1982398631705605},
  /* 34 */ {0.1965616322328226},
  /* 35 */ {0.1949590218937863},
  /* 36 */ {0.1934264036172708},
  /* 37 */ {0.1919587200065601},
  /* 38 */ {0.1905514124267734},
  /* 39 */ {0.1892003595168700},
  /* 40 */ {0.1879018247091076},
  /* 41 */ {0.1866524112389434},
  /* 42 */ {0.1854490234153689},
  /* 43 */ {0.1842888331487062},
  /* 44 */ {0.1831692509136336},
  /* 45 */ {0.1820879004699383},
  /* 46 */ {0.1810425967800402},
  /* 47 */ {0.1800313266566926},
  /* 48 */ {0.1790522317510414},
  /* 49 */ {0.1781035935540111},
  /* 50 */ {0.1771838201355579},
  /* 51 */ {0.1762914343888821},
  /* 52 */ {0.1754250635819545},
  /* 53 */ {0.1745834300480449},
  /* 54 */ {0.1737653428714400},
  /* 55 */ {0.1729696904450771},
  /* 56 */ {0.1721954337940981},
  /* 57 */ {0.1714416005739134},
  /* 58 */ {0.1707072796637201},
  /* 59 */ {0.1699916162869140},
  /* 60 */ {0.1692938075987814},
  /* 61 */ {0.1686130986895011},
  /* 62 */ {0.1679487789570419},
  /* 63 */ {0.1673001788101741},
  /* 64 */ {0.1666666666666667},
  /* 65 */ {0.1660476462159378},
  /* 66 */ {0.1654425539190583},
  /* 67 */ {0.1648508567221603},
  /* 68 */ {0.1642720499620502},
  /* 69 */ {0.1637056554452156},
  /* 70 */ {0.1631512196835108},
  /* 71 */ {0.1626083122716342},
  /* 72 */ {0.1620765243931223},
  /* 73 */ {0.1615554674429964},
  /* 74 */ {0.1610447717564444},
  /* 75 */ {0.1605440854340214},
  /* 76 */ {0.1600530732548213},
  /* 77 */ {0.1595714156699382},
  /* 78 */ {0.1590988078692941},
  /* 79 */ {0.1586349589155960},
  /* 80 */ {0.1581795909397823},
  /* 81 */ {0.1577324383928644},
  /* 82 */ {0.1572932473495469},
  /* 83 */ {0.1568617748594410},
  /* 84 */ {0.1564377883420715},
  /* 85 */ {0.1560210650222250},
  /* 86 */ {0.1556113914024939},
  /* 87 */ {0.1552085627701551},
  /* 88 */ {0.1548123827357682},
  /* 89 */ {0.1544226628011101},
  /* 90 */ {0.1540392219542636},
  /* 91 */ {0.1536618862898642},
  /* 92 */ {0.1532904886526781},
  /* 93 */ {0.1529248683028321},
  /* 94 */ {0.1525648706011593},
  /* 95 */ {0.1522103467132434},
  /* 96 */ {0.1518611533308632},
  /* 97 */ {0.1515171524096389},
  /* 98 */ {0.1511782109217764},
  /* 99 */ {0.1508442006228941},
  /* 100 */ {0.1505149978319906},
  /* 101 */ {0.1501904832236880},
  /* 102 */ {0.1498705416319474},
  /* 103 */ {0.1495550618645152},
  /* 104 */ {0.1492439365274121},
  /* 105 */ {0.1489370618588283},
  /* 106 */ {0.1486343375718350},
  /* 107 */ {0.1483356667053617},
  /* 108 */ {0.1480409554829326},
  /* 109 */ {0.1477501131786861},
  /* 110 */ {0.1474630519902391},
  /* 111 */ {0.1471796869179852},
  /* 112 */ {0.1468999356504447},
  /* 113 */ {0.1466237184553111},
  /* 114 */ {0.1463509580758620},
  /* 115 */ {0.1460815796324244},
  /* 116 */ {0.1458155105286054},
  /* 117 */ {0.1455526803620167},
  /* 118 */ {0.1452930208392429},
  /* 119 */ {0.1450364656948130},
  /* 120 */ {0.1447829506139581},
  /* 121 */ {0.1445324131589439},
  /* 122 */ {0.1442847926987864},
  /* 123 */ {0.1440400303421672},
  /* 124 */ {0.1437980688733776},
  /* 125 */ {0.1435588526911310},
  /* 126 */ {0.1433223277500932},
  /* 127 */ {0.1430884415049874},
  /* 128 */ {0.1428571428571428},
  /* 129 */ {0.1426283821033600},
  /* 130 */ {0.1424021108869747},
  /* 131 */ {0.1421782821510107},
  /* 132 */ {0.1419568500933153},
  /* 133 */ {0.1417377701235801},
  /* 134 */ {0.1415209988221527},
  /* 135 */ {0.1413064939005528},
  /* 136 */ {0.1410942141636095},
  /* 137 */ {0.1408841194731412},
  /* 138 */ {0.1406761707131039},
  /* 139 */ {0.1404703297561400},
  /* 140 */ {0.1402665594314587},
  /* 141 */ {0.1400648234939879},
  /* 142 */ {0.1398650865947379},
  /* 143 */ {0.1396673142523192},
  /* 144 */ {0.1394714728255649},
  /* 145 */ {0.1392775294872041},
  /* 146 */ {0.1390854521985406},
  /* 147 */ {0.1388952096850913},
  /* 148 */ {0.1387067714131417},
  /* 149 */ {0.1385201075671774},
  /* 150 */ {0.1383351890281539},
  /* 151 */ {0.1381519873525671},
  /* 152 */ {0.1379704747522905},
  /* 153 */ {0.1377906240751463},
  /* 154 */ {0.1376124087861776},
  /* 155 */ {0.1374358029495937},
  /* 156 */ {0.1372607812113589},
  /* 157 */ {0.1370873187823978},
  /* 158 */ {0.1369153914223921},
  /* 159 */ {0.1367449754241439},
  /* 160 */ {0.1365760475984821},
  /* 161 */ {0.1364085852596902},
  /* 162 */ {0.1362425662114337},
  /* 163 */ {0.1360779687331669},
  /* 164 */ {0.1359147715670014},
  /* 165 */ {0.1357529539050150},
  /* 166 */ {0.1355924953769864},
  /* 167 */ {0.1354333760385373},
  /* 168 */ {0.1352755763596663},
  /* 169 */ {0.1351190772136599},
  /* 170 */ {0.1349638598663645},
  /* 171 */ {0.1348099059658080},
  /* 172 */ {0.1346571975321549},
  /* 173 */ {0.1345057169479844},
  /* 174 */ {0.1343554469488779},
  /* 175 */ {0.1342063706143054},
  /* 176 */ {0.1340584713587979},
  /* 177 */ {0.1339117329233981},
  /* 178 */ {0.1337661393673756},
  /* 179 */ {0.1336216750601996},
  /* 180 */ {0.1334783246737591},
  /* 181 */ {0.1333360731748201},
  /* 182 */ {0.1331949058177136},
  /* 183 */ {0.1330548081372441},
  /* 184 */ {0.1329157659418126},
  /* 185 */ {0.1327777653067443},
  /* 186 */ {0.1326407925678156},
  /* 187 */ {0.1325048343149731},
  /* 188 */ {0.1323698773862368},
  /* 189 */ {0.1322359088617821},
  /* 190 */ {0.1321029160581950},
  /* 191 */ {0.1319708865228925},
  /* 192 */ {0.1318398080287045},
  /* 193 */ {0.1317096685686114},
  /* 194 */ {0.1315804563506306},
  /* 195 */ {0.1314521597928493},
  /* 196 */ {0.1313247675185968},
  /* 197 */ {0.1311982683517524},
  /* 198 */ {0.1310726513121843},
  /* 199 */ {0.1309479056113158},
  /* 200 */ {0.1308240206478128},
  /* 201 */ {0.1307009860033912},
  /* 202 */ {0.1305787914387386},
  /* 203 */ {0.1304574268895465},
  /* 204 */ {0.1303368824626505},
  /* 205 */ {0.1302171484322746},
  /* 206 */ {0.1300982152363760},
  /* 207 */ {0.1299800734730872},
  /* 208 */ {0.1298627138972530},
  /* 209 */ {0.1297461274170591},
  /* 210 */ {0.1296303050907487},
  /* 211 */ {0.1295152381234257},
  /* 212 */ {0.1294009178639407},
  /* 213 */ {0.1292873358018581},
  /* 214 */ {0.1291744835645007},
  /* 215 */ {0.1290623529140715},
  /* 216 */ {0.1289509357448472},
  /* 217 */ {0.1288402240804449},
  /* 218 */ {0.1287302100711566},
  /* 219 */ {0.1286208859913518},
  /* 220 */ {0.1285122442369443},
  /* 221 */ {0.1284042773229231},
  /* 222 */ {0.1282969778809442},
  /* 223 */ {0.1281903386569819},
  /* 224 */ {0.1280843525090381},
  /* 225 */ {0.1279790124049077},
  /* 226 */ {0.1278743114199984},
  /* 227 */ {0.1277702427352035},
  /* 228 */ {0.1276667996348261},
  /* 229 */ {0.1275639755045533},
  /* 230 */ {0.1274617638294791},
  /* 231 */ {0.1273601581921740},
  /* 232 */ {0.1272591522708010},
  /* 233 */ {0.1271587398372755},
  /* 234 */ {0.1270589147554692},
  /* 235 */ {0.1269596709794558},
  /* 236 */ {0.1268610025517973},
  /* 237 */ {0.1267629036018709},
  /* 238 */ {0.1266653683442337},
  /* 239 */ {0.1265683910770258},
  /* 240 */ {0.1264719661804097},
  /* 241 */ {0.1263760881150453},
  /* 242 */ {0.1262807514205999},
  /* 243 */ {0.1261859507142915},
  /* 244 */ {0.1260916806894653},
  /* 245 */ {0.1259979361142023},
  /* 246 */ {0.1259047118299582},
  /* 247 */ {0.1258120027502338},
  /* 248 */ {0.1257198038592741},
  /* 249 */ {0.1256281102107963},
  /* 250 */ {0.1255369169267456},
  /* 251 */ {0.1254462191960791},
  /* 252 */ {0.1253560122735751},
  /* 253 */ {0.1252662914786691},
  /* 254 */ {0.1251770521943144},
  /* 255 */ {0.1250882898658681},
  /* 256 */ {0.1250000000000000},
};

MPFR_COLD_FUNCTION_ATTR void
mpfr_assert_fail (const char *filename, int linenum,
                  const char *expr)
{
  if (filename != NULL && filename[0] != '\0')
    {
      fprintf (stderr, "%s:", filename);
      if (linenum != -1)
        fprintf (stderr, "%d: ", linenum);
    }
  fprintf (stderr, "MPFR assertion failed: %s\n", expr);
  abort();
}

/* Performing a concentration for these indirect functions may be
   good for performance since branch prediction for indirect calls
   is not well supported by a lot of CPU's (typically they can only
   predict a limited number of indirections). */
MPFR_HOT_FUNCTION_ATTR void *
mpfr_allocate_func (size_t alloc_size)
{
  void * (*allocate_func) (size_t);
  void * (*reallocate_func) (void *, size_t, size_t);
  void   (*free_func) (void *, size_t);

  MPFR_ASSERTD (alloc_size > 0);
  /* Always calling with the 3 arguments smooths branch prediction. */
  mp_get_memory_functions (&allocate_func, &reallocate_func, &free_func);
  return (*allocate_func) (alloc_size);
}

MPFR_HOT_FUNCTION_ATTR void *
mpfr_reallocate_func (void * ptr, size_t old_size, size_t new_size)
{
  void * (*allocate_func) (size_t);
  void * (*reallocate_func) (void *, size_t, size_t);
  void   (*free_func) (void *, size_t);

  MPFR_ASSERTD (new_size > 0);
  /* Always calling with the 3 arguments smooths branch prediction. */
  mp_get_memory_functions (&allocate_func, &reallocate_func, &free_func);
  return (*reallocate_func) (ptr, old_size, new_size);
}

MPFR_HOT_FUNCTION_ATTR void
mpfr_free_func (void *ptr, size_t size)
{
  void * (*allocate_func) (size_t);
  void * (*reallocate_func) (void *, size_t, size_t);
  void   (*free_func) (void *, size_t);

  MPFR_ASSERTD (size > 0);
  /* Always calling with the 3 arguments smooths branch prediction. */
  mp_get_memory_functions (&allocate_func, &reallocate_func, &free_func);
  (*free_func) (ptr, size);
}

void *
mpfr_tmp_allocate (struct tmp_marker **tmp_marker, size_t size)
{
  struct tmp_marker *head;

  head = (struct tmp_marker *)
    mpfr_allocate_func (sizeof (struct tmp_marker));
  head->ptr = mpfr_allocate_func (size);
  head->size = size;
  head->next = *tmp_marker;
  *tmp_marker = head;
  return head->ptr;
}

void
mpfr_tmp_free (struct tmp_marker *tmp_marker)
{
  struct tmp_marker *t;

  while (tmp_marker != NULL)
    {
      t = tmp_marker;
      mpfr_free_func (t->ptr, t->size);
      tmp_marker = t->next;
      mpfr_free_func (t, sizeof (struct tmp_marker));
    }
}

#endif /* Have gmp-impl.h */
