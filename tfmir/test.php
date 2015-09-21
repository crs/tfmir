<?php

  echo getcwd();
  echo get_current_user();
  $old_umask = umask(0);
  mkdir('./uploads/allaccess', 0775);
  //chown('./uploads/allaccess', get_current_user());
  //chmod('./uploads/allaccess', 0775);
?>
