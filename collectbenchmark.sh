/usr/bin/bash
echo -n name; cat * | head -n1 ; ls | parallel 'echo -n {}\\t ;tail -n+2 {}'


