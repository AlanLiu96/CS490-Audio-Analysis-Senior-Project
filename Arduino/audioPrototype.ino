//Audio out with 38.5kHz sampling rate and interrupts
//by Amanda Ghassaei
//https://www.instructables.com/id/Arduino-Audio-Input/
//Sept 2012

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
*/

int incomingAudio;
//int t = 0;

void setup(){

  cli();//disable interrupts
  
  //set up continuous sampling of analog pin 0
  
  //clear ADCSRA and ADCSRB registers
  ADCSRA = 0;
  ADCSRB = 0;
  
  ADMUX |= (1 << REFS0); //set reference voltage
  ADMUX |= (1 << ADLAR); //left align the ADC value- so we can read highest 8 bits from ADCH register only
  
  ADCSRA |= (1 << ADPS2) | (1 << ADPS0); //set ADC clock with 32 prescaler- 16mHz/32=500kHz
  ADCSRA |= (1 << ADATE); //enabble auto trigger
  ADCSRA |= (1 << ADIE); //enable interrupts when measurement complete
  ADCSRA |= (1 << ADEN); //enable ADC
  ADCSRA |= (1 << ADSC); //start ADC measurements

    //set timer0 interrupt at 40kHz
  TCCR0A = 0;// set entire TCCR0A register to 0
  TCCR0B = 0;// same for TCCR0B
  TCNT0  = 0;//initialize counter value to 0
  // set compare match register for 40khz increments
  OCR0A = 49;// = (16*10^6) / (40000*8) - 1 (must be <256)
  // turn on CTC mode
  TCCR0A |= (1 << WGM01);
  // Set CS11 bit for 8 prescaler
  TCCR0B |= (1 << CS11); 
  // enable timer compare interrupt
  TIMSK0 |= (1 << OCIE0A);
  
  sei();//enable interrupts

  //if you want to add other things to setup(), do it here
  Serial.begin(9600);

  //set digital pins 0-7 as outputs
  for (int i=0;i<8;i++){
    pinMode(i,OUTPUT);
  }
  
}

ISR(TIMER0_COMPA_vect){ //40kHz interrupt routine
  PORTD = incomingAudio;//send audio wave to DAC, centered around (127/255)*5 = 2.5V
}

ISR(ADC_vect) {//when new ADC value ready
  incomingAudio = ADCH;//update the variable incomingAudio with new value from A0 (between 0 and 255)
}

void loop(){
  Serial.println(incomingAudio);
  //do other stuff here
}


