How to contribute
-----------------

We are happy to accept pull requests.

We use a modified gnu coding style. The code can and should be formatted
using the indent tool like this: ::

       indent -gnu -fc1 -i4 -bli0 -nut -cdb -sc -bap -l80 *.c && rm -rf *~

Block comments should look like this and precede functions etc.: ::

      /*
       * I am a block comment!
       */

For one-lined comments, use :code:`//`
