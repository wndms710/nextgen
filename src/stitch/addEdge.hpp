{
  
  const int ied_new = new_edge();
  nooed[ied_new][0] = current_ino;
  nooed[ied_new][1] = next_ino;
  
  switch(current_id) {
    // =========================================
    // z-edges...
    // =========================================
  case 0:
    faoed[ied_new][0] = -2;
    faoed[ied_new][1] = -4;
    break;
  case 1:
    faoed[ied_new][0] = -4;
    faoed[ied_new][1] = -3;
    break;
  case 2:
    faoed[ied_new][0] = -3;
    faoed[ied_new][1] = -5;
    break;
  case 3:
    faoed[ied_new][0] = -5;
    faoed[ied_new][1] = -2;
    break;
    // =========================================
    // x-edges...
    // =========================================
  case 4:
    faoed[ied_new][0] = -4;
    faoed[ied_new][1] = -6;
    break;
  case 5:
    faoed[ied_new][0] = -6;
    faoed[ied_new][1] = -5;
    break;
  case 6:
    faoed[ied_new][0] = -5;
    faoed[ied_new][1] = -7;
    break;
  case 7:
    faoed[ied_new][0] = -7;
    faoed[ied_new][1] = -4;
    break;
    // =========================================
    // y-edges...
    // =========================================
  case 8:
    faoed[ied_new][0] = -6;
    faoed[ied_new][1] = -2;
    break;
  case 9:
    faoed[ied_new][0] = -2;
    faoed[ied_new][1] = -7;
    break;
  case 10:
    faoed[ied_new][0] = -7;
    faoed[ied_new][1] = -3;
    break;
  case 11:
    faoed[ied_new][0] = -3;
    faoed[ied_new][1] = -6;
    break;
  default:
    assert(0);
  }

}
	    
