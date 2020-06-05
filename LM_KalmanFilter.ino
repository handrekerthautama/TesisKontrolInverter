String dataPC[4], dat;
int data[4], j=0;
float e[2], esum[] = {5.39, 0}, tmp[2], v[2], dz[]={25, 4865, 0, 0}, y[2], state[4], tmp_state[]={0.25, 48.65, 0, 0};

void setup() {

  Serial.begin(9600);
  
}

void loop() {
 if(Serial.available() > 0){
    for(int i=0; i<2; i=i+1){
      dataPC[i] = Serial.readStringUntil(',');
      tmp[i] = dataPC[i].toInt();
      y[i] = tmp[i] / 1000;
    }

    dz[0] = dz[0] + 0.0559 * (y[0] - tmp_state[0]);
    dz[1] = dz[1] - 0.5000 * (y[0] - tmp_state[0]);
    dz[2] = dz[2] + 0.0559 * (y[1] - tmp_state[2]);
    dz[3] = dz[3] - 0.5000 * (y[1] - tmp_state[2]);

    if(j == 100){
      state[0] = 225.5;
    }
    
    for(int i=0; i<4; i=i+1){
      state[i] = state[i] + dz[i] * 0.01;
    }
    
    for(int i=0; i<2; i=i+1){
       e[i] = ((1-i) * 230 - state[i*2]);
       esum[i] = esum[i] + e[i] * 0.01;
       v[i] = 2236.1 * esum[i] - 365.3 * state[i*2] - 27 * state[i*2+1];
    }

    dz[0] = state[1];
    dz[1] = v[0];
    dz[2] = state[3];
    dz[3] = v[1];

    for(int i=0; i<4; i=i+1){
      tmp_state[i] = state[i] + dz[i] * 0.01;
    }
    
    dat = String(v[0]) + ',' + String(v[1]) + ',' + String(tmp_state[1]) + ',' + String(tmp_state[3]);
    Serial.print(dat);
    j = j+1;
   }
 }
