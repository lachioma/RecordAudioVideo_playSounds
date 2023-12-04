
const int pinLed1 =  8; // Led1 for Sounds
const int pinLed2 =  9; // Led2 for Visuals
const int pinLed3 = 10; // Led that copies both Led1 and Led2
const int pinLed4 = 13; // Built-in/debug Led that copies both Led1 and Led2

char receivedChar;
boolean newData = false;

const int durLed1      = 1000; // how many msec keep the LED1 on upon trigger
const int durLed2      = 2000; // how many msec keep the LED2 on upon trigger
const int delayLed1    =    0; // how many msec wait after trigger before switching LED1 on
const int delayLed2    =    0; // how many msec wait after trigger before switching LED2 on
boolean LED1isON = false;
boolean LED2isON = false;
unsigned long t0 = millis(); // will store time LED was turned on
unsigned long tF = millis(); // will store time TTL_frame was turned on

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
  while (!Serial) {
    ; //
  }
  // Confirm connection
  Serial.println("Arduino online");
  pinMode(pinLed1, OUTPUT);
  digitalWrite(pinLed1, LOW);
  pinMode(pinLed2, OUTPUT);
  digitalWrite(pinLed2, LOW);
  pinMode(pinLed3, OUTPUT);
  digitalWrite(pinLed3, LOW);
  pinMode(pinLed4, OUTPUT);
  digitalWrite(pinLed4, LOW);
}

void loop() {
  // put your main code here, to run repeatedly:
  newData = false;
  recvOneChar();

  if (newData == true && receivedChar == 'S') {
    t0 = millis();
    while ((millis()-t0) < delayLed1) {
      ; // do nothing, just wait
    }
    digitalWrite(pinLed1, HIGH);
    //digitalWrite(pinLed2, HIGH);
    digitalWrite(pinLed3, HIGH);
    digitalWrite(pinLed4, HIGH);
    t0 = millis();
    //Serial.println("LED turned on");
    newData = false;
    LED1isON = true;
  }
  else if (newData == true && receivedChar == 'V') {
    t0 = millis();
    while ((millis()-t0) < delayLed2) {
      ; // do nothing, just wait
    }
    //digitalWrite(pinLed1, HIGH);
    digitalWrite(pinLed2, HIGH);
    digitalWrite(pinLed3, HIGH);
    digitalWrite(pinLed4, HIGH);
    t0 = millis();
    //Serial.println("LED turned on");
    newData = false;
    LED2isON = true;
  }
  else if (newData == true && receivedChar == '0') {
    digitalWrite(pinLed1, LOW);
    digitalWrite(pinLed2, LOW);
    digitalWrite(pinLed3, LOW);
    digitalWrite(pinLed4, LOW);
    //Serial.println("LED turned off");
    newData = false;
    LED1isON = false;
    LED2isON = false;
  }

  else if (newData == true && 
          (receivedChar != '0' && receivedChar != 'S' && receivedChar != 'V') ) {
    digitalWrite(pinLed1, LOW);
    digitalWrite(pinLed2, LOW);
    digitalWrite(pinLed3, LOW);
    digitalWrite(pinLed4, LOW);
    //Serial.println("LED turned off");
    newData = false;
    LED1isON = false;
    LED2isON = false;
  }

  if (LED1isON == true && ((millis()-t0) >= durLed1)) {
    digitalWrite(pinLed1, LOW);
    //digitalWrite(pinLed2, LOW);
    digitalWrite(pinLed3, LOW);
    digitalWrite(pinLed4, LOW);
    LED1isON = false;
  }
  if (LED2isON == true && ((millis()-t0) >= durLed2)) {
    //digitalWrite(pinLed1, LOW);
    digitalWrite(pinLed2, LOW);
    digitalWrite(pinLed3, LOW);
    digitalWrite(pinLed4, LOW);
    LED1isON = false;
  }
}

void recvOneChar() {
    // When the LED and TTL are Off, just wait here until a message is received.
    // Wait here only when the LED and TTL are Off, otherwise you will not be able to switch them Off.
    // Waiting here probably saves just a few microseconds.
    while (Serial.available() < 1 && LED1isON == false && LED2isON == false) {
      ;
      // delayMicroseconds(1);
    }
    if (Serial.available() > 0) {
        receivedChar = Serial.read();
        newData = true;
    }
}