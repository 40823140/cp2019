ä¸­æ–‡// ¤U¦C Dart µ{¦¡, §Q¥Î Runge Kutta ­¡¥N¹Bºâªk, ¸Ñ±`·L¤À¤èµ{¦¡
// ³] t ¬°®É¶¡, x «h³]¬°ª«Åéªº¦ì²¾
// °²³]­n¸Ñ F=ma ªº³æ¤@½è¶q¥[¤W¼uÂ® (±`¼Æ¬° k) »PÂHº¢ªý¥§ (±`¼Æ¬° b)
// f ¬°ªu¦ì²¾¤è¦Vªº¬I¤O
// dx/dt = v, dv/dt = (f-kx-bv)/m
// dx / dt = (t - x)/2, °_©l­È t0=0, x0=1, ¨D t=2 ®Éªº x ­È
//
// ¤wª¾°_©l­È t0 »P x0 «á, ¥i¥H§Q¥Î¤U¦C rungeKutta ¨ç¦¡, ¥H
// h ¬°¨C¨B¶¥¼W¶q­È, ¨D dxdt ±`·L¤À¤èµ{¦¡¥ô¤@ t ªº¹ïÀ³­È x
// ©w¸q¨ç¦¡ rungeKutta, ¦@¦³¥|­Ó¿é¤JÅÜ¼Æ
// ª«Åé½è¶q
const num m = 1;
// ¹ï½è¶qªº¬I¤O f
const num f = 0.0;
// ¼uÂ®«Y¼Æ
const num k = 1;
// ªý¥§«Y¼Æ
const num b = 1;

// ©I¥s¹Bºâ®É, »Ý­n°_©l®É¶¡, ²×ÂI®É¶¡, ¦ì²¾°_©l­È»P³t«×°_©l­È, ¼W¶q h
rungeKutta(t0, x0, v0, t, h) {
  // §Q¥Î¨B¶¥¼W¶q­È h »P t ªº°_©l¤Î²×ÂI­È
  // ­pºâ»Ý­n­¡¥Nªº¦¸¼Æ n
  int n = ((t - t0) / h).toInt();
  // «Å§i x ¬°Âù¯BÂI¼Æ, ¥B³]¬°°_©l­È x0
  double x = x0;
  // «Å§i v ¬°Âù¯BÂI¼Æ, ¥B³]¬°°_©l­È v0
  double v = v0;

  // ¼ÒÀÀ¹Bºâ«e, ¦C¥X°_©l±ø¥ó
  // ¥u¦C¨ì¤p¼ÆÂI²Ä¤T¦ì, ®É¶¡¡B¦ì²¾»P³t«×¥H \t  ¹j¶}, \t ¥Nªí´¡¤J tab ²Å¸¹, ¥i±N¸ê®Æ½Æ»s¨ì Excel ¶i¦æÃ¸¹Ï
  print("${t0.toStringAsFixed(3)} \t ${x.toStringAsFixed(3)} \t ${v.toStringAsFixed(3)}");

  // §Q¥Î¤wª¾ªº t0, x0, t ²×ÂI­È»P¨B¶¥¼W¶q­È h, ­¡¥N¨D x ¹ïÀ³­È
  // ¯Á¤Þ­È i ±N¨C¦¸¼W¶q 1, ±q i=1 °õ¦æ for Àô°é¦Ü i=n
  for (int i = 1; i <= n; i++) {
    // ±N¦¹¶¥¬qªº t »P x ­È¥N¤J dxdt »P dvdt ¨ç¦¡¨D¤U¦C¥|­Ó¯BÂIÅÜ¼Æ­È
    // ¦]¬°¥²¶·¨â­Ó¨ç¦¡½¢¦X, ¥²¶·¦P®É­pºâ
    double xk1 = h * dxdt(t0, x, v);
    double vk1 = h * dvdt(t0, x, v);
    double xk2 = h * dxdt(t0 + 0.5 * h, x + 0.5 * xk1, v + 0.5 * vk1);
    double vk2 = h * dvdt(t0 + 0.5 * h, x + 0.5 * xk1, v + 0.5 * vk1);
    double xk3 = h * dxdt(t0 + 0.5 * h, x + 0.5 * xk2, v + 0.5 * vk2);
    double vk3 = h * dvdt(t0 + 0.5 * h, x + 0.5 * xk2, v + 0.5 * vk2);
    double xk4 = h * dxdt(t0 + h, x + xk3, 