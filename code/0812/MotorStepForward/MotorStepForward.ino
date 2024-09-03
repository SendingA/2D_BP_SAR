#define DIR 10
#define PLS 11
#define DIRM 8
#define PLSM 9

#define DIR2 6
#define PLS2 7
#define DIRM2 4
#define PLSM2 5

void setup() {
  pinMode(PLS, OUTPUT);
  pinMode(DIR, OUTPUT);
  pinMode(PLSM, OUTPUT);
  pinMode(DIRM, OUTPUT);

  pinMode(PLS2, OUTPUT);
  pinMode(DIR2, OUTPUT);
  pinMode(PLSM2, OUTPUT);
  pinMode(DIRM2, OUTPUT);

  digitalWrite(PLS, HIGH);
  digitalWrite(DIR, HIGH);
  digitalWrite(PLSM, LOW);
  digitalWrite(DIRM, LOW);

  digitalWrite(PLS2, HIGH);
  digitalWrite(DIR2, HIGH);
  digitalWrite(PLSM2, LOW);
  digitalWrite(DIRM2, LOW);
}

// DIR X轴，HIGH 右正，LOW 左负
// DIR2 Y轴，HIGH 下正，LOW 上负


void loop() {
  
  delay(1000);

//  PLS, DIR, 水平电机, HIGH远离电机，LOW靠近电机
//  PLS2, DID2, 竖直电机，LOW向下, HIGH向上
//  数值5对应1cm
//  1cm对应0.5s运动时间


//   scanStep(PLS2, DIR2, HIGH, 25);

//  scanStep(PLS, DIR, LOW, 200);

  for(int i_h = 0; i_h < 30; i_h++){
    scanStep(PLS, DIR, HIGH, 5);
    delay(1500);
  }

//  delay(2000);
//
//  for(int i_v = 0; i_v < 30; i_v ++){
//    scanStep(PLS2, DIR2, HIGH, 5);
//    delay(1000);
//  }
//  
  while(1);

}

void scanStep(int plsPin, int dirPin, int dirHL, int bigSteps) {
  digitalWrite(dirPin, dirHL);
  int steps = 1600/2;

  for(int s=0; s<bigSteps; s++){
    for(int x=0; x<steps; x++){
      digitalWrite(plsPin, HIGH);
      delayMicroseconds(625/10);
      digitalWrite(plsPin, LOW);
      delayMicroseconds(625/10);
      }
  }
}
