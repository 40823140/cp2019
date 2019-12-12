中文// �U�C Dart �{��, �Q�� Runge Kutta ���N�B��k, �ѱ`�L����{��
// �] t ���ɶ�, x �h�]�����骺�첾
// ���]�n�� F=ma ����@��q�[�W�u® (�`�Ƭ� k) �P�H������ (�`�Ƭ� b)
// f ���u�첾��V���I�O
// dx/dt = v, dv/dt = (f-kx-bv)/m
// dx / dt = (t - x)/2, �_�l�� t0=0, x0=1, �D t=2 �ɪ� x ��
//
// �w���_�l�� t0 �P x0 ��, �i�H�Q�ΤU�C rungeKutta �禡, �H
// h ���C�B���W�q��, �D dxdt �`�L����{�����@ t �������� x
// �w�q�禡 rungeKutta, �@���|�ӿ�J�ܼ�
// �����q
const num m = 1;
// ���q���I�O f
const num f = 0.0;
// �u®�Y��
const num k = 1;
// �����Y��
const num b = 1;

// �I�s�B���, �ݭn�_�l�ɶ�, ���I�ɶ�, �첾�_�l�ȻP�t�װ_�l��, �W�q h
rungeKutta(t0, x0, v0, t, h) {
  // �Q�ΨB���W�q�� h �P t ���_�l�β��I��
  // �p��ݭn���N������ n
  int n = ((t - t0) / h).toInt();
  // �ŧi x �����B�I��, �B�]���_�l�� x0
  double x = x0;
  // �ŧi v �����B�I��, �B�]���_�l�� v0
  double v = v0;

  // �����B��e, �C�X�_�l����
  // �u�C��p���I�ĤT��, �ɶ��B�첾�P�t�ץH \t  �j�}, \t �N���J tab �Ÿ�, �i�N��ƽƻs�� Excel �i��ø��
  print("${t0.toStringAsFixed(3)} \t ${x.toStringAsFixed(3)} \t ${v.toStringAsFixed(3)}");

  // �Q�Τw���� t0, x0, t ���I�ȻP�B���W�q�� h, ���N�D x ������
  // ���ޭ� i �N�C���W�q 1, �q i=1 ���� for ����� i=n
  for (int i = 1; i <= n; i++) {
    // �N�����q�� t �P x �ȥN�J dxdt �P dvdt �禡�D�U�C�|�ӯB�I�ܼƭ�
    // �]��������Ө禡���X, �����P�ɭp��
    double xk1 = h * dxdt(t0, x, v);
    double vk1 = h * dvdt(t0, x, v);
    double xk2 = h * dxdt(t0 + 0.5 * h, x + 0.5 * xk1, v + 0.5 * vk1);
    double vk2 = h * dvdt(t0 + 0.5 * h, x + 0.5 * xk1, v + 0.5 * vk1);
    double xk3 = h * dxdt(t0 + 0.5 * h, x + 0.5 * xk2, v + 0.5 * vk2);
    double vk3 = h * dvdt(t0 + 0.5 * h, x + 0.5 * xk2, v + 0.5 * vk2);
    double xk4 = h * dxdt(t0 + h, x + xk3, 