import fit

specials = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
    'Saturn', 'Uranus', 'Neptune']
src_data = {
    #             RA          DEC         FLUX  FREQ, INDEX ANGSIZE
    #'Moon'  : (None        , None,           0, .150,  0.10 ),
    #'Mercury': (None        , None,          0, .150,  0.10 ),
    #'Venus':   (None        , None,          0, .150,  0.10 ),
    #'Mars':    (None        , None,          0, .150,  0.10 ),
    #'Jupiter': (None        , None,          0, .150,  0.10 ),
    #'Saturn' : (None        , None,          0, .150,  0.10 ),
    #'Uranus' : (None        , None,          0, .150,  0.10 ),
    #'Neptune': (None        , None,          0, .150,  0.10 ),
    'Sun'  : (None, None, 57000, .150,  2.00, 4.6e-3),
    'for':   ('03:22:41.7', '-37:12:30',   170, .150, -0.8 , 0.),
    'pic'  : ('05:19:49.7', '-45:46:45',   452, .160, -0.8 , 0.),
    'hyd'  : ('09:18:05.7', '-12:05:44',  1860, .160, -2.30, 0.),
    'cen'  : ('13:25:27.6', '-43:01:09', 7104, .160, -0.52, 3.5e-3),
    'her' : ('16:51:08.15', '4:59:33.3', 300.0, 0.159, -1, 0.000669),
    'sgr'  : ('17:45:40.0', '-29:00:28',  121, .160, -4.21, 0.),
    # Fluxes from Miriad:
    'crab' : ('05:34:32.0', '+22:00:52',  1838, .150, -0.30, 0.),
    'vir'  : ('12:30:49.4', '+12:23:28',   1446, .150, -0.86, 0.),
    'cyg'  : ('19:59:28.3', '+40:44:02', 10900, .150, -0.69, (0.,0.,0.)),
    'cas'  : ('23:23:27.94', '+58:48:42.4',  9160.0, 0.150, -0.73, 0.000000),
    # 3C + 3CR Catalogs (position corrected with NED)
    '3c002' : ('00:06:22.59', '-00:04:24.7', 16.5, 0.159, -0.84, 0.000727),
    '3c009' : ('00:20:25.22', '+15:40:54.6', 15.5, 0.159, -0.29, 0.000436),
    '3c010' : ('00:25:13.99', '+64:08:39.5', 110.0, 0.159, 1.75, 0.000960),
    '3c013' : ('00:34:14.55', '+39:24:16.7', 12.0, 0.159, -1.18, 0.002909),
    '3c014' : ('00:36:06.50', '+18:37:58.6', 13.5, 0.159, -3.11, 0.002909),
    '3c015' : ('00:37:04.15', '-01:09:08.1', 21.5, 0.159, -3.49, 0.005818),
    '3c016' : ('00:37:44.57', '+13:19:55.0', 10.5, 0.159, 0.81, 0.000145),
    '3c017' : ('00:38:20.51', '-02:07:40.7', 21.5, 0.159, -0.21, 0.005818),
    '3c018' : ('00:40:50.47', '+10:03:22.7', 16.5, 0.159, 2.94, 0.005818),
    '3c019' : ('00:40:55.04', '+33:10:08.0', 13.0, 0.159, -3.26, 0.001745),
    '3c020' : ('00:43:08.84', '+52:03:33.8', 30.0, 0.159, 2.77, 0.000145),
    '3c022' : ('00:50:56.30', '+51:12:03.0', 13.5, 0.159, 0.63, 0.000291),
    '3c027' : ('00:56:01.06', '+68:22:30.4', 22.0, 0.159, 0.77, 0.000291),
    '3c028' : ('00:55:50.33', '+26:24:34.4', 19.5, 0.159, -2.94, 0.000436),
    '3c029' : ('00:57:34.91', '-01:23:27.9', 17.0, 0.159, -1.11, 0.000436),
    '3c031' : ('01:07:24.96', '+32:24:45.2', 10.5, 0.159, 3.45, 0.000364),
    '3c033' : ('01:08:52.86', '+13:20:13.8', 58.0, 0.159, -1.49, 0.000291),
    '3c034' : ('01:10:18.67', '+31:47.20.5', 12.5, 0.159, -1.13, 0.000364),
    '3c035' : ('01:12:02.23', '+49:28:35.2',  9.0, 0.159, 2.91, 0.000465),
    '3c036' : ('01:17:59.48', '+45:36:21.8',  8.0, 0.159, 1.52, 0.000582),
    '3c040' : ('01:25:59.83', '-01:20:34.4', 26.0, 0.159, -0.71, 0.000727),
    '3c041' : ('01:26:44.39', '+33:13:11.2',  8.5, 0.159, 1.44, 0.000436),
    '3c042' : ('01:28:30.33', '+29:02:59.4', 12.0, 0.159, -0.38, 0.000364),
    '3c043' : ('01:29:59.81', '+23:38:20.3', 16.5, 0.159, -3.20, 0.000436),
    '3c044' : ('01:31:21.76', '+06:23:40.8',  8.0, 0.159, 2.82, 0.000291),
    '3c046' : ('01:35:28.32', '+37:54:05.6', 12.5, 0.159, -2.43, 0.000364),
    '3c047' : ('01:36:24.40', '+20:57:27.0', 27.0, 0.159, -2.66, 0.000233),
    '3c048' : ('01:37:41.30', '+33:09:35.1', 50.0, 0.159, -0.55, 0.000145),
    '3c049' : ('01:41:09.16', '+13:53:28.0', 10.0, 0.159, 1.24, 0.000582),
    '3c052' : ('01:48:28.98', '+53:32:35.4',  9.5, 0.159, 0.45, 0.000436),
    '3c054' : ('01:55:30.16', '+43:45:55.4', 11.5, 0.159, -1.24, 0.001745),
    '3c055' : ('01:57:10.51', '+28:51:37.6',  7.5, 0.159, 9.12, 0.000364),
    '3c058' : ('02:05:35.24', '+64:49:34.8', 13.0, 0.159, 5.43, 0.000654),
    '3c063' : ('02:20:54.25', '-01:56:51.8', 26.0, 0.159, -4.87, 0.005818),
    '3c065' : ('02:23:43.19', '+40:00:52.5', 18.5, 0.159, 3.35, 0.000538),
    '3c066' : ('02:22:25.55', '+43:00:47.7', 28.0, 0.159, 1.46, 0.000480),
    '3c067' : ('02:24:12.30', '+27:50:11.5', 10.0, 0.159, -0.00, 0.000436),
    '3c069' : ('02:38:02.35', '+59:11:50.0', 23.0, 0.159, -0.00, 0.001745),
    '3c071' : ('02:42:40.71', '-00:00:47.8', 11.0, 0.159, 1.81, 0.005818),
    '3c075' : ('02:57:41.57', '+06:01:28.8', 38.0, 0.159, -4.45, 0.000524),
    '3c078' : ('03:08:26.22', '+04:06:39.3', 17.0, 0.159, -1.11, 0.000291),
    '3c079' : ('03:10:00.09', '+17:05:58.3', 34.0, 0.159, -3.09, 0.000218),
    '3c084' : ('03:19:48.16', '+41:30:42.1', 50.0, 0.159, 1.31, 0.000335),
    '3c086' : ('03:27:19.36', '+55:20:28.1', 19.0, 0.159, 0.89, 0.000364),
    '3c088' : ('03:27:54.19', '+02:33:42.0',  9.0, 0.159, 5.10, 0.000465),
    '3c089' : ('03:34:15.57', '-01:10:56.4', 19.5, 0.159, -0.47, 0.000364),
    '3c091' : ('03:37:43.35', '+50:45:52.8', 14.0, 0.159, -1.00, 0.002182),
    '3c093' : ('03:43:30.01', '+04:57:48.6', 11.5, 0.159, -1.69, 0.000436),
    '3c098' : ('03:58:54.43', '+10:26:03.0', 41.0, 0.159, -0.00, 0.000364),
    '3c099' : ('04:01:07.63', '+00:36:32.9', 14.5, 0.159, -3.29, 0.000436),
    '3c103' : ('04:08:03.34', '+43:00:24.4', 29.0, 0.159, -1.31, 0.005818),
    '3c105' : ('04:07:16.47', '+03:42:25.8', 12.5, 0.159, 1.62, 0.000509),
    '3c107' : ('04:12:22.63', '-00:59:32.4', 11.5, 0.159, -0.39, 0.001745),
    '3c109' : ('04:13:40.37', '+11:12:13.8', 19.5, 0.159, -0.00, 0.000218),
    '3c111' : ('04:18:21.28', '+38:01:35.8', 60.0, 0.159, -0.45, 0.000422),
    '3c114' : ('04:20:22.17', '+17:53:55.2', 12.5, 0.159, -1.98, 0.000436),
    '3c119' : ('04:32:36.50', '+41:38:28.4', 14.5, 0.159, 0.30, 0.000218),
    '3c123' : ('04:37:04.37', '+29:40:13.8', 204.0, 0.159, -1.36, 0.005818),
    '3c124' : ('04:41:59.11', '+01:21:01.9',  8.5, 0.159, 0.99, 0.000582),
    '3c125' : ('04:46:17.86', '+39:45:03.0', 12.5, 0.159, -2.43, 0.000436),
    '3c129' : ('04:49:09.07', '+45:00:39.2', 21.5, 0.159, -0.21, 0.000393),
    '3c130' : ('04:52:52.84', '+52:04:47.1',  9.5, 0.159, 2.43, 0.000436),
    '3c131' : ('04:53:23.33', '+31:29:24.3', 16.0, 0.159, -1.51, 0.000291),
    '3c132' : ('04:56:43.05', '+22:49:22.2', 16.5, 0.159, -2.46, 0.000145),
    '3c133' : ('05:02:58.55', '+25:16:24.7', 23.0, 0.159, -1.46, 0.006545),
    '3c134' : ('05:04:42.19', '+38:06:11.4', 85.0, 0.159, -2.24, 0.005818),
    '3c135' : ('05:14:08.36', '+00:56:32.5', 12.0, 0.159, 2.55, 0.000538),
    '3c137' : ('05:19:32.42', '+50:54:32.1',  8.5, 0.159, 0.51, 0.005818),
    '3c138' : ('05:21:09.88', '+16:38:22.1', 19.5, 0.159, -0.47, 0.000436),
    '3c141' : ('05:26:44.20', '+32:50:23.0', 20.5, 0.159, -4.04, 0.005818),
    #'3c144' : ('05:34:31.97', '+22:00:52.1', 1500.0, 0.159, -0.49, 0.000727),
    '3c147' : ('05:42:36.14', '+49:51:07.2', 63.0, 0.159, -0.73, 0.000291),
    '3c152' : ('06:04:28.63', '+20:21:21.7', 12.5, 0.159, -1.13, 0.001745),
    '3c153' : ('06:09:32.50', '+48:04:15.4', 15.0, 0.159, 1.37, 0.005818),
    '3c154' : ('06:13:50.14', '+26:04:36.7', 26.0, 0.159, -2.78, 0.000145),
    '3c157' : ('06:16:37.41', '+22:31:54.0', 270.0, 0.159, -2.23, 0.004218),
    '3c158' : ('06:21:41.06', '+14:32:11.2', 21.5, 0.159, -3.80, 0.001745),
    '3c165' : ('06:43:07.40', '+23:19:02.6', 12.5, 0.159, -1.13, 0.000364),
    '3c166' : ('06:45:24.10', '+21:21:51.2', 16.0, 0.159, -1.84, 0.000364),
    '3c171' : ('06:55:14.81', '+54:09:00.1', 30.0, 0.159, -2.35, 0.006545),
    '3c172' : ('07:02:08.07', '+25:13:46.3', 17.0, 0.159, -1.72, 0.000364),
    '3c173' : ('07:02:17.60', '+37:57:19.5', 15.5, 0.159, -3.88, 0.000436),
    '3c175' : ('07:13:02.40', '+11:46:14.7', 23.5, 0.159, -3.41, 0.006545),
    '3c177' : ('07:25:09.00', '+15:36:00.0', 12.5, 0.159, -2.43, 0.000436),
    '3c180' : ('07:27:04.50', '-02:04:42.3', 15.0, 0.159, -0.61, 0.000436),
    '3c181' : ('07:28:10.30', '+14:37:36.2', 14.0, 0.159, -0.66, 0.000364),
    '3c184' : ('07:39:24.47', '+70:23:10.9', 17.0, 0.159, -3.86, 0.000364),
    '3c186' : ('07:44:17.45', '+37:53:17.1', 14.0, 0.159, -0.32, 0.001745),
    '3c187' : ('07:45:04.47', '+02:00:08.1',  9.0, 0.159, 1.78, 0.000436),
    '3c190' : ('08:01:33.51', '+14:14:42.4', 12.0, 0.159, -0.77, 0.000582),
    '3c191' : ('08:04:47.93', '+10:15:23.2', 11.0, 0.159, -0.41, 0.005818),
    '3c192' : ('08:05:35.00', '+24:09:50.0', 17.0, 0.159, 1.22, 0.000364),
    '3c194' : ('08:10:03.62', '+42:28:04.3', 12.5, 0.159, -2.91, 0.000436),
    '3c196' : ('08:13:36.03', '+48:13:02.6', 66.0, 0.159, -0.99, 0.001745),
    '3c198' : ('08:22:31.95', '+05:57:06.8', 16.0, 0.159, 0.27, 0.000538),
    '3c200' : ('08:27:25.38', '+29:18:45.4', 10.5, 0.159, 1.89, 0.000436),
    '3c204' : ('08:37:44.96', '+65:13:34.9',  9.0, 0.159, 0.48, 0.000436),
    '3c205' : ('08:39:06.46', '+57:54:17.1', 11.0, 0.159, 1.13, 0.005818),
    '3c207' : ('08:40:47.59', '+13:12:23.6', 10.0, 0.159, -0.00, 0.000436),
    '3c208' : ('08:53:08.85', '+13:52:55.5', 24.5, 0.159, -3.77, 0.000364),
    '3c210' : ('08:58:09.96', '+27:50:51.6', 10.0, 0.159, 0.43, 0.001745),
    '3c212' : ('08:58:41.50', '+14:09:44.0', 21.5, 0.159, -2.90, 0.002182),
    '3c215' : ('09:06:31.90', '+16:46:11.4', 15.5, 0.159, -3.88, 0.000436),
    '3c216' : ('09:09:33.50', '+42:53:46.5', 23.5, 0.159, -2.12, 0.001745),
    '3c217' : ('09:08:50.58', '+37:48:19.2', 10.5, 0.159, 1.18, 0.000364),
    '3c219' : ('09:21:08.63', '+45:38:57.4', 42.0, 0.159, 0.41, 0.000145),
    '3c222' : ('09:36:32.02', '+04:22:10.3', 12.0, 0.159, -2.55, 0.000582),
    '3c223' : ('09:39:52.74', '+35:53:58.2', 14.0, 0.159, 0.31, 0.000364),
    '3c225' : ('09:42:12.00', '+13:48:52.0', 19.5, 0.159, 2.20, 0.000436),
    '3c226' : ('09:44:16.37', '+09:46:19.2', 11.0, 0.159, 1.81, 0.000364),
    '3c227' : ('09:47:45.12', '+07:25:20.6', 50.0, 0.159, -5.14, 0.000407),
    '3c228' : ('09:50:10.74', '+14:19:58.1', 17.0, 0.159, -0.54, 0.005818),
    '3c230' : ('09:51:58.82', '-00:01:27.2', 31.0, 0.159, -3.45, 0.001745),
    '3c231' : ('09:55:52.72', '+69:40:45.8', 12.0, 0.159, 0.71, 0.000218),
    '3c234' : ('10:01:49.55', '+28:47:09.3', 30.0, 0.159, -0.30, 0.005818),
    '3c236' : ('10:06:01.74', '+34:54:10.4', 12.0, 0.159, -2.55, 0.000436),
    '3c237' : ('10:08:00.03', '+07:30:16.3', 21.5, 0.159, 0.97, 0.001745),
    '3c238' : ('10:11:00.38', '+06:24:39.7', 18.5, 0.159, -3.13, 0.000364),
    '3c239' : ('10:11:45.41', '+46:28:19.8', 15.0, 0.159, -1.98, 0.000364),
    '3c241' : ('10:21:54.55', '+21:59:30.1', 13.0, 0.159, -2.32, 0.006545),
    '3c245' : ('10:42:44.60', '+12:03:31.3', 12.0, 0.159, -2.07, 0.005818),
    '3c247' : ('10:58:58.78', '+43:01:23.1',  9.5, 0.159, 4.34, 0.000218),
    '3c249' : ('11:02:03.85', '-01:16:17.4', 14.5, 0.159, 2.85, 0.000364),
    '3c250' : ('11:08:52.13', '+25:00:54.6', 14.0, 0.159, -2.14, 0.000436),
    '3c252' : ('11:11:33.08', '+35:40:41.6', 15.0, 0.159, -1.98, 0.000291),
    '3c254' : ('11:14:38.71', '+40:37:20.3', 21.5, 0.159, -1.10, 0.001745),
    '3c255' : ('11:19:25.24', '-03:02:51.5', 15.0, 0.159, 0.84, 0.000364),
    '3c256' : ('11:20:43.05', '+23:27:54.9', 11.5, 0.159, -1.69, 0.001745),
    '3c257' : ('11:23:09.17', '+05:30:19.5', 11.0, 0.159, -1.78, 0.000582),
    '3c258' : ('11:24:43.53', '+19:19:11.6', 10.0, 0.159, -0.93, 0.000582),
    '3c263' : ('11:39:57.04', '+65:47:49.4', 10.0, 0.159, 2.32, 0.002182),
    '3c264' : ('11:45:05.01', '+19:36:22.7', 37.0, 0.159, -3.83, 0.000393),
    '3c265' : ('11:45:28.99', '+31:33:49.4', 30.0, 0.159, -4.53, 0.000145),
    '3c266' : ('11:45:43.37', '+49:46:08.2', 14.0, 0.159, -3.91, 0.000436),
    '3c267' : ('11:49:56.52', '+12:47:18.9', 14.5, 0.159, -0.63, 0.000509),
    '3c270' : ('12:19:23.22', '+05:49:30.8', 20.0, 0.159, 6.98, 0.000684),
    '3c272' : ('12:24:28.52', '+42:06:36.3', 10.5, 0.159, -0.89, 0.000436),
    '3c273' : ('12:29:06:70', '+02:03:08.6', 79.0, 0.159, -1.46, 0.001745),
    #'3c274' : ('12:30:49.42', '+12:23:28.0', 1100.0, 0.159, -1.11, 0.000684),
    '3c275' : ('12:42:19.66', '-04:46:20.7', 18.0, 0.159, -5.66, 0.000436),
    '3c277' : ('12:51:43.58', '+50:34:24.9', 12.0, 0.159, 0.36, 0.000582),
    '3c280' : ('12:56:57.09', '+47:20:19.6', 25.0, 0.159, -1.98, 0.005818),
    '3c284' : ('13:11:04.67', '+27:28:07.7', 10.0, 0.159, 0.43, 0.000436),
    '3c285' : ('13:21:17.82', '+42:35:15.2', 11.0, 0.159, -0.41, 0.000436),
    '3c286' : ('13:31:08.29', '+30:30:33.0', 30.0, 0.159, -3.16, 0.000218),
    '3c287' : ('13:30:37.69', '+25:09:10.9', 29.0, 0.159, -6.45, 0.001745),
    '3c288' : ('13:38:49.87', '+38:51:09.2', 15.0, 0.159, -0.30, 0.005818),
    '3c289' : ('13:45:26.37', '+49:46:32.6', 12.0, 0.159, -2.07, 0.000364),
    '3c293' : ('13:52:17.84', '+31:26:46.5', 12.5, 0.159, -0.36, 0.000364),
    '3c294' : ('14:06:44.02', '+34:11:25.1', 12.5, 0.159, -2.91, 0.005818),
    '3c295' : ('14:11:20.65', '+52:12:09.0', 74.0, 0.159, -0.12, 0.001745),
    '3c296' : ('14:16:52.95', '+10:48:26.6', 10.0, 0.159, 1.98, 0.000436),
    '3c297' : ('14:17:24.00', '-04:00:47.5', 14.5, 0.159, -3.75, 0.000000),
    '3c298' : ('14:19:08.18', '+06:28:34.8', 61.0, 0.159, -2.89, 0.000582),
    '3c299' : ('14:21:05.58', '+41:44:48.5', 10.5, 0.159, 0.41, 0.000291),
    '3c300' : ('14:23:01.03', '+19:35:17.4', 15.0, 0.159, 0.57, 0.000291),
    '3c303' : ('14:43:02.76', '+52:01:37.2',  9.0, 0.159, 2.91, 0.005818),
    '3c305' : ('14:49:21.57', '+63:16:14.0', 15.0, 0.159, -0.93, 0.001745),
    '3c310' : ('15:04:57.12', '+26:00:58.5', 72.0, 0.159, -3.05, 0.000233),
    '3c313' : ('15:11:00.01', '+07:51:50.0', 21.0, 0.159, 2.55, 0.000364),
    '3c315' : ('15:13:40.07', '+26:07:31.2', 26.0, 0.159, -3.51, 0.000291),
    '3c317' : ('15:16:44.49', '+07:01:16.6', 55.0, 0.159, -2.18, 0.005818),
    '3c318' : ('15:20:05.45', '+20:16:05.8', 14.5, 0.159, -3.29, 0.002182),
    '3c319' : ('15:24:05.50', '+54:28:14.6', 16.5, 0.159, 0.26, 0.000364),
    '3c320' : ('15:31:25.37', '+35:33:40.0',  8.0, 0.159, 3.21, 0.005818),
    '3c321' : ('15:31:43.45', '+24:04:19.1', 15.0, 0.159, -1.62, 0.005818),
    '3c322' : ('15:35:01.23', '+55:36:52.9', 11.5, 0.159, -1.69, 0.000436),
    '3c323' : ('15:41:45.53', '+60:15:35.1',  9.0, 0.159, 0.48, 0.000393),
    '3c324' : ('15:49:48.90', '+21:25:38.1', 18.0, 0.159, -3.97, 0.005818),
    '3c325' : ('15:49:58.57', '+62:41:20.7', 15.0, 0.159, -3.16, 0.005818),
    '3c326' : ('15:52:09.15', '+20:05:23.7', 12.5, 0.159, -1.98, 0.000247),
    '3c327' : ('16:02:27.35', '+01:57:56.2', 34.0, 0.159, 1.44, 0.000538),
    '3c330' : ('16:09:36.61', '+65:56:43.6', 24.0, 0.159, -0.00, 0.000291),
    '3c332' : ('16:17:42.52', '+32:22:34.8', 16.0, 0.159, -3.73, 0.000436),
    '3c334' : ('16:20:21.92', '+17:36:24.0', 16.0, 0.159, -4.16, 0.005818),
    '3c336' : ('16:24:39.09', '+23:45:12.2', 13.5, 0.159, -0.00, 0.005818),
    '3c337' : ('16:28:52.84', '+44:19:05.1',  8.5, 0.159, 6.14, 0.000436),
    '3c338' : ('16:28:38.48', '+39:33:05.6', 49.0, 0.159, -1.58, 0.000291),
    '3c340' : ('16:29:36.93', '+23:20:14.4',  9.5, 0.159, 1.30, 0.000436),
    '3c341' : ('16:28:04.05', '+27:41:43.0',  8.5, 0.159, 1.87, 0.000436),
    '3c343' : ('16:34:33.79', '+62:45:35.8', 18.0, 0.159, -5.66, 0.000436),
    '3c345' : ('16:42:58.81', '+39:48:37.0',  9.0, 0.159, 0.93, 0.001745),
    '3c346' : ('16:43:48.60', '+17:15:49.5', 15.5, 0.159, -3.04, 0.002182),
    #'3c348' : ('16:51:08.15', '+04:59:33.3', 300.0, 0.159, 0.71, 0.000465),
    '3c349' : ('16:59:29.54', '+47:02:44.1', 13.5, 0.159, -1.04, 0.005818),
    '3c351' : ('17:04:41.37', '+60:44:30.5', 15.0, 0.159, -2.75, 0.000436),
    '3c352' : ('17:10:44.10', '+46:01:28.6', 12.0, 0.159, -0.77, 0.001745),
    '3c353' : ('17:20:28.15', '-00:58:46.8', 180.0, 0.159, 1.07, 0.000582),
    '3c356' : ('17:24:19.02', '+50:57:40.3', 14.0, 0.159, 0.31, 0.000218),
    '3c357' : ('17:28:18.47', '+31:46:14.4',  9.0, 0.159, -0.00, 0.000436),
    '3c368' : ('18:05:06.36', '+11:01:32.5', 13.5, 0.159, 0.32, 0.000509),
    '3c371' : ('18:06:50.68', '+69:49:28.1',  9.0, 0.159, 0.48, 0.000436),
    '3c380' : ('18:29:31.78', '+48:44:46.2', 70.0, 0.159, -1.82, 0.002909),
    '3c381' : ('18:33:46.29', '+47:27:02.7', 14.5, 0.159, -1.31, 0.002909),
    '3c382' : ('18:35:03.39', '+32:41:46.8', 18.0, 0.159, 0.93, 0.000262),
    '3c386' : ('18:38:12.85', '+17:11:49.3', 27.0, 0.159, -1.04, 0.000262),
    '3c388' : ('18:44:02.40', '+45:33:29.7', 22.5, 0.159, -0.20, 0.005818),
    '3c389' : ('18:46:18.63', '-03:19:44.5', 17.0, 0.159, 1.87, 0.000480),
    '3c390' : ('18:45:37.62', '+09:53:44.7', 22.5, 0.159, -1.50, 0.005818),
    '3c391' : ('18:49:21.57', '-00:55:32.3', 27.0, 0.159, -1.04, 0.000291),
    '3c392' : ('18:56:10.00', '+01:19:00.0', 680.0, 0.159, -4.70, 0.002327),
    '3c394' : ('18:59:23.36', '+12:59:12.1', 18.5, 0.159, -3.83, 0.000145),
    '3c396' : ('19:04:04.48', '+05:27:12.4', 24.5, 0.159, -0.95, 0.000436),
    '3c397' : ('19:07:32.90', '+07:08:33.0', 29.0, 0.159, -0.31, 0.000436),
    '3c398' : ('19:11:05.33', '+09:05:38.9', 43.0, 0.159, 2.18, 0.000480),
    '3c400' : ('19:22:58.00', '+14:11:50.3', 25.0, 0.159, 27.22, 0.008727),
    '3c401' : ('19:40:25.12', '+60:41:34.9', 22.0, 0.159, -1.30, 0.005818),
    '3c402' : ('19:41:46.00', '+50:35:44.9', 15.0, 0.159, -1.27, 0.000393),
    '3c403' : ('19:52:15.80', '+02:30:24.5', 23.0, 0.159, 1.09, 0.000218),
    #'3c405' : ('19:59:28.35', '+40:44:02.1', 8600.0, 0.159, -0.53, 0.000000),
    '3c409' : ('20:14:27.59', '+23:34:52.9', 102.0, 0.159, -2.96, 0.004363),
    '3c410' : ('20:20:06.56', '+29:42:14.2', 36.0, 0.159, -1.04, 0.001745),
    '3c411' : ('20:22:08.45', '+10:01:11.7', 16.0, 0.159, -1.18, 0.004363),
    '3c418' : ('20:38:37.03', '+51:19:12.7', 16.0, 0.159, -2.19, 0.000436),
    '3c424' : ('20:48:12.03', '+07:01:17.5', 16.0, 0.159, -2.19, 0.005818),
    '3c428' : ('21:08:22.39', '+49:36:37.6', 19.5, 0.159, -0.23, 0.000364),
    '3c430' : ('21:18:19.11', '+60:48:07.4', 100.0, 0.159, -10.97, 0.000262),
    '3c431' : ('21:18:52.46', '+49:36:59.1', 31.0, 0.159, -3.45, 0.000218),
    '3c432' : ('21:22:46.25', '+17:04:38.3', 13.5, 0.159, -1.42, 0.005818),
    '3c433' : ('21:23:44.53', '+25:04:11.9', 62.0, 0.159, -1.56, 0.004363),
    '3c434' : ('21:23:16.35', '+15:48:06.1', 10.5, 0.159, -0.00, 0.000509),
    '3c435' : ('21:29:06.10', '+07:33:06.4', 12.5, 0.159, -2.43, 0.005818),
    '3c436' : ('21:44:11.73', '+28:10:18.4', 21.0, 0.159, -2.98, 0.005818),
    '3c437' : ('21:47:25.10', '+15:20:37.5', 16.0, 0.159, -3.73, 0.005818),
    '3c438' : ('21:55:52.29', '+38:00:29.6', 43.0, 0.159, -1.33, 0.005818),
    '3c441' : ('22:06:04.91', '+29:29:20.0', 12.5, 0.159, 0.35, 0.005818),
    '3c442' : ('22:14:46.90', '+13:50:24.0', 33.0, 0.159, -4.44, 0.000000),
    '3c445' : ('22:23:49.57', '-02:06:12.4', 20.5, 0.159, 1.02, 0.000276),
    '3c449' : ('22:31:20.90', '+39:21:48.0', 11.5, 0.159, 1.42, 0.000436),
    '3c452' : ('22:45:48.77', '+39:41:15.7', 50.0, 0.159, -0.18, 0.000465),
    '3c454' : ('22:51:34.74', '+18:48:40.1', 11.5, 0.159, -2.17, 0.000436),
    '3c455' : ('22:55:03.83', '+13:13:34.2', 15.0, 0.159, -1.27, 0.000364),
    '3c456' : ('23:12:28.07', '+09:19:29.2', 10.0, 0.159, 2.32, 0.000436),
    '3c458' : ('23:12:54.41', '+05:16:46.0', 12.5, 0.159, -1.98, 0.000320),
    '3c459' : ('23:16:35.23', '+04:05:18.1', 25.0, 0.159, -1.13, 0.000436),
    '3c460' : ('23:21:28.53', '+23:46:48.4', 13.0, 0.159, -2.78, 0.000364),
    #'3c461' : ('23:23:27.94', '+58:48:42.4', 13000.0, 0.159, -1.48, 0.000000),
    '3c465' : ('23:38:29.38', '+27:01:53.2', 50.0, 0.159, -3.16, 0.000582),
    '3c470' : ('23:58:35.34', '+44:04:38.9',  8.0, 0.159, 1.04, 0.000364),
}

def get_src(s, fixedbody=fit.RadioFixedBody, special=fit.RadioSpecial):
    """Return a source created out of the parameters in the dictionary srcs.
    Can pass your own RadioFixedBody or RadioSpecial subclasses to use."""
    if not type(s) == str: return s
    ra, dec, st, mfreq, index, srcshape = src_data[s]
    if s in specials:
        return special(s, st, mfreq=mfreq, 
            index=index, srcshape=srcshape)
    else:
        return fixedbody(ra, dec, janskies=st, mfreq=mfreq, 
            index=index, name=s, srcshape=srcshape)

def get_catalog(srcs=None, cutoff=None, 
        fixedbody=fit.RadioFixedBody, special=fit.RadioSpecial):
    """Return a source catalog created out of the parameters in the 
    dictionary srcs.  Can pass your own RadioFixedBody or RadioSpecial 
    subclasses to use."""
    if srcs is None:
        if cutoff is None: srcs = src_data.keys()
        else: srcs = [s for s in src_data if src_data[s][2] > cutoff]
    srcs = [get_src(s, fixedbody=fixedbody, special=special) for s in srcs]
    return fit.SrcCatalog(srcs)
